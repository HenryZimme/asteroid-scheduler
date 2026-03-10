# asteroid-scheduler

**Phase-gap-aware multi-site observing window scheduler for asteroid rotational lightcurve campaigns.**

Given a known or candidate rotation period, this tool identifies when and where uncovered phase intervals (gaps) in a lightcurve are observable across a set of ground-based sites. It queries live ephemerides from the Minor Planet Center, computes per-night visibility windows using dynamically calculated astronomical twilight, and evaluates lunar sky brightness via the Krisciunas-Schaefer (1991) physical model. Output is printed to the console and saved to CSV.

Developed as part of an ongoing photometric campaign on asteroid **7605 Cindygraber** at Phillips Academy Observatory (MPC I12).

---

## Motivation

Single-site lightcurve campaigns on slowly-rotating asteroids frequently produce phase gaps — intervals of the rotation period that are perpetually unobservable from a fixed longitude due to the relationship between synodic drift and the local night window. Closing a phase gap requires coordinated observations from geographically separated sites during the specific UTC intervals when the gap phase is visible.

No existing tool in the Python small-body astronomy ecosystem addresses this problem directly.

---

## Features

- **MPC ephemeris query** via `astroquery.mpc.MPC.get_ephemeris` — live RA/Dec, magnitude, distance, and altitude for any numbered asteroid at any site and step size
- **Phase gap scheduling** — given T0 (BJD_TDB+LTTC), period, and gap phase boundaries, computes both daily gap occurrences and their observable windows at each site
- **General visibility mode** — works without lightcurve parameters as a standalone multi-site observing window calculator
- **Astronomical twilight** — computed dynamically per site and date using pyephem (Sun at −18°); no hardcoded values
- **Lunar sky brightness** — full Krisciunas & Schaefer (1991) PASP 103, 1033 physical model; reports V-band sky brightness (mag/arcsec²) at window start and midpoint with qualitative quality labels
- **Automatic timezone detection** — IANA timezone strings resolved via `timezonefinder`; DST handled automatically with no seasonal edits
- **BJD_TDB + light travel time correction** — phase calculation uses barycentric time with light travel time subtracted per Eastman et al. (2010)
- **CSV output** — separate files for the full ephemeris table and per-site/per-night observing windows with all moon metrics
- **iCalendar export** — optional `.ics` files for import into calendar apps

---

## Requirements

```
astroquery     >= 0.4.6
astropy        >= 5.0
ephem          >= 4.1
timezonefinder >= 6.0
tzdata
```

Install with:

```bash
pip install astroquery astropy ephem timezonefinder tzdata
```

The script is designed for use in Google Colab but runs in any standard Python 3.9+ environment.

---

## Usage

All user-facing configuration lives in **Cell 1** of `asteroid_scheduler.py`. No other cell needs to be modified for a new target or set of sites.

### Minimal configuration (ephemeris only)

```python
ASTEROID_ID   = '7605'
ASTEROID_NAME = 'Cindygraber'
DATE_START    = datetime(2026, 3, 1)
N_DAYS        = 10
EPHEM_STEP    = '1h'
USE_LIGHTCURVE = False

SITES = {
    'G40 Canary 1': {
        'lat': 28.29970, 'lon': -16.5082, 'elev': 2390,
        'mpc_code': 'G40',
    },
}
PRIMARY_SITE = 'G40 Canary 1'
```

Timezone is resolved automatically from latitude/longitude. Running with `USE_LIGHTCURVE = False` produces a nightly visibility window for each site (twilight, altitude > `MIN_ALT_DEG`) with lunar sky brightness, without any phase filtering.

### Full phase gap configuration

```python
USE_LIGHTCURVE = True
T0_JD         = 2461084.655574   # BJD_TDB+LTTC from lightcurve fit
PERIOD_H      = 11.939            # rotation period in hours
GAP_START_PH  = 0.72              # phase where your coverage ends
GAP_END_PH    = 1.00              # phase 1.0 = 0.0
```

`T0_JD` must be in Barycentric Julian Date (BJD_TDB) with the one-way light travel time correction already applied. If your epoch is in JD_UTC, convert it first using the LTTC calculator included in Cell 1 or via [Eastman et al. (2010)](https://astroutils.astronomy.osu.edu/time/bjdconvert.html).

The scheduler computes two gap occurrences per night (separated by one rotation period), intersects each with the dark/visible window at every site, and reports observable windows, altitudes, and sky conditions for Gap 1 and Gap 2 independently.

### Adding sites

```python
SITES = {
    'G40 Canary 1': {
        'lat': 28.29970, 'lon': -16.5082, 'elev': 2390,
        'mpc_code': 'G40',
    },
    'E62 Australia 1': {
        'lat': -31.2816709, 'lon': 149.0801825, 'elev': 805,
        'mpc_code': 'E62',
    },
    'W88 Chile 2': {
        'lat': -33.269, 'lon': -70.53, 'elev': 1492,
        'mpc_code': 'W88',
    },
    'I12 Phillips Academy Observatory': {
        'lat': 42.647611, 'lon': -71.129, 'elev': 80,
        'mpc_code': 'I12',
    },
}
```

Any MPC observatory code is valid. Timezone offsets and DST are resolved automatically from coordinates. Add or remove entries freely — the scheduler loops over all sites.

---

## Output

### Console

```
============================================================================================
  ASTEROID OBSERVING SCHEDULER
  Target   : 7605 Cindygraber
  Dates    : 2026-03-01  to  2026-03-10  (10 nights)
  Mode     : phase gap filter ON
  T0       : JD 2461084.655574   Period: 11.939 h
  Gap      : phase 0.72 – 1.0  (200.6 min = 3h 20m)  drift 7.32 min/day earlier
  Alt min  : > 20°   step: 1 min
  Sky model: Krisciunas-Schaefer (1991)   k = 0.172
============================================================================================

  SITE: G40 CANARY 1
  Location  : 28.2997° N   16.5082° W   2390 m
  MPC code  : G40
  Time zone : UTC+0
  Dark sky  : 21.9 mag/arcsec²

  -- GAP 2 -----------------------------------------------------------------------

  Date     Gap window         Min    UTC window             Local window              Altitude       Limit
  ------------------------------------------------------------------------------------------
  03-01   gap 22:19-01:39    201    22:19 – 01:39 UTC    22:19 – 01:39 local (UTC+0)    64° → 54°     gap ends
           Moon illum  97%  start:  17.15 mag/"  Poor     sep   12°  moon   62° up  mid:  17.16 mag/"  Poor     sep   12°  moon   77° up
```

### CSV files

**`ephem_<ID>.csv`** — full MPC ephemeris at the queried step resolution. Columns include Date, RA, Dec, V, Delta, r, Altitude, Azimuth, and all columns returned by MPC.

**`windows_<ID>.csv`** — one row per site/date/gap with columns:

| Column | Description |
|--------|-------------|
| `site` | Site name |
| `date` | UTC date |
| `gap_label` | GAP 1 or GAP 2 |
| `gap_window` | Phase gap UTC window |
| `obs_start_utc` / `obs_end_utc` | Observable window in UTC |
| `obs_start_local` / `obs_end_local` | Observable window in local time |
| `time_offset` | UTC offset string |
| `duration_min` | Observable duration in minutes |
| `alt_start` / `alt_end` | Asteroid altitude at window start/end |
| `limit` | Constraining factor (gap ends / dawn / altitude) |
| `start_v_sky` | Sky brightness at window start (mag/arcsec²) |
| `start_sky_label` | Excellent / Good / Fair / Poor |
| `start_illum_pct` | Moon illumination % at start |
| `start_moon_alt` | Moon altitude at start |
| `start_separation` | Moon-asteroid separation at start |
| `mid_*` | Same five columns evaluated at window midpoint |

---

## Sky brightness model

Lunar sky brightness is computed using the physical model of [Krisciunas & Schaefer (1991)](https://ui.adsabs.harvard.edu/abs/1991PASP..103.1033K), which accounts for:

- Moon phase angle derived from illumination fraction
- Moon V-band magnitude as a function of phase angle
- Scattering function f(ρ) for moon-target angular separation
- Atmospheric extinction at both moon and target airmasses (Rozenberg 1966)
- Site-specific dark sky baseline

The result is the total V-band sky surface brightness in mag/arcsec². Quality labels use configurable thresholds (defaults: Excellent ≥ 21.5, Good ≥ 20.5, Fair ≥ 19.5, Poor < 19.5).

---

## Validation: 7605 Cindygraber

The tool was developed and validated against an active photometric campaign on **7605 Cindygraber** (main belt, H = 13.6) observed from Phillips Academy Observatory (MPC I12) beginning late 2025.

| Parameter | Value |
|-----------|-------|
| Rotation period | 11.939 h |
| T0 (BJD_TDB+LTTC) | 2461084.655574 |
| Phase gap | 0.72 – 1.00 (200.6 min) |
| Gap drift | 7.32 min/day earlier |
| Campaign window | 2026 March 1–10 |

**Site coverage:**

| Site | Code | Gap | Observable nights |
|------|------|-----|-------------------|
| Canary 1 | G40 | Gap 2 (22:19 UTC) | All 10, full 201 min |
| Australia 1 | E62 | Gap 1 (10:22 UTC) | All 10, 168–201 min |

Moon interference (illumination 55–100%, separation 11–22° at campaign start) degrades sky brightness to 17–18 mag/arcsec² for March 1–5. Conditions improve to Excellent (≥ 21.5) from March 6 as the moon sets before the observing windows open.

---

## Configuration reference

| Parameter | Type | Description |
|-----------|------|-------------|
| `ASTEROID_ID` | str | MPC number or packed designation |
| `ASTEROID_NAME` | str | Display name only |
| `DATE_START` | datetime | First night of campaign (UTC midnight) |
| `N_DAYS` | int | Number of nights |
| `EPHEM_STEP` | str | MPC query step: `'1d'`, `'6h'`, `'1h'`, `'30min'` |
| `ALT_STEP_MINUTES` | int | Altitude stepping resolution for window calculation |
| `SITES` | dict | Site configs: lat, lon, elev, mpc_code |
| `PRIMARY_SITE` | str | Site used for MPC query |
| `MIN_ALT_DEG` | float | Minimum asteroid altitude threshold |
| `DARK_SKY_MAG` | float or dict | Site dark sky V brightness (mag/arcsec²) |
| `EXTINCTION_K` | float or dict | V-band extinction coefficient (mag/airmass) |
| `SKY_EXCELLENT` / `SKY_GOOD` / `SKY_FAIR` | float | Quality label thresholds |
| `USE_LIGHTCURVE` | bool | Enable phase gap filtering |
| `T0_JD` | float | Epoch of phase 0.0 (BJD_TDB+LTTC) |
| `PERIOD_H` | float | Rotation period in hours |
| `GAP_START_PH` | float | Phase where existing coverage ends |
| `GAP_END_PH` | float | Phase where gap closes (typically 1.0) |
| `SAVE_CSV` | bool | Write output CSV files |
| `CSV_EPHEM` / `CSV_WINDOWS` | str | Output filenames |

---

## Cell structure

The script is organized into numbered cells for use in Google Colab. Only **Cell 1** (user configuration) should ever be edited.

| Cell | Contents |
|------|----------|
| 1 | User configuration — edit for each new asteroid/campaign |
| 2 | Imports |
| 3 | MPC ephemeris query |
| 4 | Site setup, altitude, and twilight utilities |
| 5 | RA/Dec interpolation from MPC table |
| 6 | BJD_TDB + light travel time correction |
| 7 | Krisciunas-Schaefer sky brightness model |
| 8 | Window calculator |
| 9 | Display and CSV output |
| 10 | Main |

---

## Citation

If you use this code in published work, please cite the software record:

> Zimmerman, H. (2026). *asteroid-scheduler: Phase-gap-aware multi-site observing window scheduler for asteroid rotational lightcurve campaigns*. Astrophysics Source Code Library. [ascl:XXXX.XXX]

A BibTeX entry is provided in `CITATION.cff`.

The Cindygraber campaign results are described in:

> Zimmerman, H. (in prep). Rotation period and taxonomy of 7605 Cindygraber. *Minor Planet Bulletin*.

---

## License

MIT License — see [LICENSE](LICENSE).

## Author

Henry Zimmerman  
Phillips Academy, Andover, MA  
MPC Observatory Code: I12
