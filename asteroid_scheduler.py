"""
asteroid ephemeris + observing window scheduler
with Krisciunas-Schaefer (1991) lunar sky brightness model
==========================================================
queries MPC via astroquery.mpc.MPC.get_ephemeris, computes visibility
windows at each site (twilight, altitude), optionally constrains to a
lightcurve phase gap, and reports lunar sky brightness at window start
and midpoint. local times shown alongside UTC using per-site offsets.

install: pip install astroquery astropy ephem
"""

# ==============================================================================
# CELL 1 — USER CONFIGURATION (only cell that changes between asteroids)
# ==============================================================================

# --- target ---
ASTEROID_ID   = '7605'
ASTEROID_NAME = 'Cindygraber'

# --- date range ---
from datetime import datetime, timedelta
DATE_START = datetime(2026, 3, 1)
N_DAYS     = 10

# --- ephemeris sampling ---
# step passed to MPC.get_ephemeris. examples: '1d', '6h', '1h', '30min'
# '1h' recommended — finer steps give better RA/Dec interpolation.
# console table shows one row per day; full data goes to CSV.
EPHEM_STEP = '1h'

# resolution for altitude/moon stepping in the window calculator (minutes)
ALT_STEP_MINUTES = 1

# --- observing sites ---
# time_offset: hours from UTC (e.g. -5 = EST, +11 = AEDT)
# update seasonally for sites that observe daylight saving time. current settings are set to March 1, 2026
SITES = {
    'G40 Canary 1': {
        'lat':         28.29970,
        'lon':        -16.5082,
        'elev':        2390,
        'mpc_code':   'G40',
        'time_offset': 0,       # WET (winter); change to +1 for summer
    },
    'W88 Chile 2': {
        'lat':         -33.269,
        'lon':         -70.53,
        'elev':        1492,
        'mpc_code':   'W88',
        'time_offset': -3,      # CLST (summer); change to -4 for winter
    },
    'E62 Australia 1': {
        'lat':         -31.2816709,
        'lon':         149.0801825,
        'elev':        805,
        'mpc_code':   'E62',
        'time_offset': 11,      # AEDT (summer); change to 10 for winter
    },
    'I12 Phillips Academy Observatory': {
        'lat':         42.647611,
        'lon':        -71.129,
        'elev':        80,
        'mpc_code':   'I12',
        'time_offset': -5,      # EST (winter); change to -4 for summer
    },
}

# site used for the MPC ephemeris query
PRIMARY_SITE = 'G40 Canary 1'

# --- altitude threshold ---
MIN_ALT_DEG = 20.0

# --- moon / sky brightness (Krisciunas-Schaefer 1991) ---
# dark sky V-band brightness at each site (mag/arcsec^2).
# typical: excellent site ~22.0, good ~21.5, average ~21.0
# either a single float (all sites) or a dict keyed by site name.
DARK_SKY_MAG = {
    'G40 Canary 1':                    21.9,
    'E62 Australia 1':                 21.7,
    'W88 Chile 2':                     21.8,
    'I12 Phillips Academy Observatory': 20.5,
}

# V-band extinction coefficient (mag/airmass)
# typical sea-level ~0.20, high-altitude good site ~0.15-0.17
EXTINCTION_K = 0.172

# sky brightness quality thresholds (mag/arcsec^2, higher = darker = better)
SKY_EXCELLENT = 21.5
SKY_GOOD      = 20.5
SKY_FAIR      = 19.5
# below SKY_FAIR → 'Poor'

# --- lightcurve parameters ---
# USE_LIGHTCURVE = False  →  full visibility window per night (no phase filter)
# USE_LIGHTCURVE = True   →  windows further restricted to the phase gap;
#                            Gap 1 and Gap 2 shown separately
USE_LIGHTCURVE = True

T0_JD        = 2461084.655574   # JDo(LTC) from lightcurve plot: epoch of phase 0.0
PERIOD_H     = 11.939            # rotation period in hours
GAP_START_PH = 0.72              # phase where your coverage ends
GAP_END_PH   = 1.00              # phase 1.0 = 0.0 (full rotation)

# --- output files ---
SAVE_CSV    = True
CSV_EPHEM   = f'ephem_{ASTEROID_ID}.csv'
CSV_WINDOWS = f'windows_{ASTEROID_ID}.csv'

# ==============================================================================
# CELL 2 — IMPORTS
# ==============================================================================

import ephem
import math
import csv
from astroquery.mpc import MPC

# ==============================================================================
# CELL 3 — EPHEMERIS QUERY
# ==============================================================================

def query_ephemeris(asteroid_id, site_code, start_date, n_days, step):
    """
    query MPC for ephemeris. returns astropy table.
    step examples: '1d', '6h', '1h', '30min'
    """
    step_lower = step.lower()
    if step_lower.endswith('d'):
        step_h = float(step_lower[:-1]) * 24
    elif step_lower.endswith('h'):
        step_h = float(step_lower[:-1])
    elif step_lower.endswith('min'):
        step_h = float(step_lower[:-3]) / 60
    else:
        raise ValueError(
            f"unrecognised step: '{step}'. use e.g. '1d', '6h', '30min'"
        )

    n_steps = int((n_days * 24) / step_h) + 1
    print(f"  querying MPC: {asteroid_id}  site {site_code}  "
          f"{start_date.strftime('%Y-%m-%d')}  {n_steps} rows x {step} ...")

    eph = MPC.get_ephemeris(
        target   = asteroid_id,
        location = site_code,
        start    = start_date.strftime('%Y-%m-%d'),
        number   = n_steps,
        step     = step,
    )
    print(f"  done — {len(eph)} rows retrieved")
    return eph

# ==============================================================================
# CELL 4 — SITE / ALTITUDE / TWILIGHT UTILITIES
# ==============================================================================

def make_obs(cfg):
    obs = ephem.Observer()
    obs.lat       = str(cfg['lat'])
    obs.lon       = str(cfg['lon'])
    obs.elevation = cfg['elev']
    obs.horizon   = '-18'   # astronomical twilight: sun at -18 deg
    obs.pressure  = 0
    return obs

def get_alt_ephem(obs, dt, ra_str, dec_str):
    a = ephem.FixedBody()
    a._ra    = ephem.hours(ra_str)
    a._dec   = ephem.degrees(dec_str)
    a._epoch = ephem.J2000
    obs.date = dt.strftime('%Y/%m/%d %H:%M:%S')
    a.compute(obs)
    return math.degrees(float(a.alt))

def evening_twilight_end(obs, date):
    """utc time when sun crosses -18 deg going down on the evening of 'date'"""
    sun = ephem.Sun()
    obs.date = date.strftime('%Y/%m/%d 12:00:00')
    return ephem.Date(obs.next_setting(sun, use_center=True)).datetime()

def morning_dawn(obs, date):
    """utc time when sun crosses -18 deg going up on the morning of 'date'"""
    sun = ephem.Sun()
    obs.date = date.strftime('%Y/%m/%d 00:01:00')
    return ephem.Date(obs.next_rising(sun, use_center=True)).datetime()

def dark_window(obs, date, gap_s=None):
    """
    return (dark_start, dark_end) = astronomical twilight end → dawn.
    uses gap_s to determine which night's dark window to use:
      gap in utc morning (hour < 12) → previous evening's twilight → this morning's dawn
      gap in utc evening/night        → this evening's twilight → next morning's dawn
    when gap_s is None (no lightcurve), defaults to evening.
    """
    if gap_s is not None and gap_s.hour < 12:
        dark_s = evening_twilight_end(obs, date - timedelta(days=1))
        dark_e = morning_dawn(obs, date)
    else:
        dark_s = evening_twilight_end(obs, date)
        dark_e = morning_dawn(obs, date + timedelta(days=1))
    return dark_s, dark_e

def utc_to_local(dt, time_offset):
    """convert utc datetime to local time using a fixed hour offset"""
    return dt + timedelta(hours=time_offset)

def offset_str(time_offset):
    """format time offset as UTC+X or UTC-X"""
    sign = '+' if time_offset >= 0 else ''
    return f"UTC{sign}{time_offset}"

# ==============================================================================
# CELL 5 — RA/Dec INTERPOLATION FROM MPC TABLE
# ==============================================================================

def build_radec_lookup(eph_table):
    """build sorted list of (datetime, ra_deg, dec_deg) from MPC astropy table"""
    lookup = []
    for row in eph_table:
        dt  = row['Date'].to_datetime()
        ra  = float(row['RA'])
        dec = float(row['Dec'])
        lookup.append((dt, ra, dec))
    return sorted(lookup, key=lambda x: x[0])

def interpolate_radec(lookup, dt):
    """
    linear interpolation of RA/Dec for a given datetime.
    returns (ra_str, dec_str) in ephem-compatible h:m:s / ±d:m:s format.
    """
    if dt <= lookup[0][0]:
        ra_deg, dec_deg = lookup[0][1], lookup[0][2]
    elif dt >= lookup[-1][0]:
        ra_deg, dec_deg = lookup[-1][1], lookup[-1][2]
    else:
        for i in range(len(lookup) - 1):
            t0, r0, d0 = lookup[i]
            t1, r1, d1 = lookup[i + 1]
            if t0 <= dt <= t1:
                frac    = (dt - t0).total_seconds() / (t1 - t0).total_seconds()
                ra_deg  = r0 + frac * (r1 - r0)
                dec_deg = d0 + frac * (d1 - d0)
                break

    ra_h   = ra_deg / 15.0
    ra_str = (f"{int(ra_h)}:{int((ra_h % 1) * 60)}:"
              f"{((ra_h % 1) * 60 % 1) * 60:.2f}")
    sign    = '+' if dec_deg >= 0 else '-'
    dec_abs = abs(dec_deg)
    dec_str = (f"{sign}{int(dec_abs)}:{int((dec_abs % 1) * 60)}:"
               f"{((dec_abs % 1) * 60 % 1) * 60:.2f}")
    return ra_str, dec_str

# ==============================================================================
# CELL 6 — KRISCIUNAS-SCHAEFER (1991) LUNAR SKY BRIGHTNESS MODEL
# ==============================================================================

def _airmass(alt_deg):
    """rozenberg (1966) airmass. returns None if more than 1 deg below horizon."""
    if alt_deg < -1.0:
        return None
    z_rad = math.radians(max(0.0, 90.0 - alt_deg))
    return 1.0 / (math.cos(z_rad) + 0.025 * math.exp(-11.0 * math.cos(z_rad)))

def _sky_label(v_sky):
    if v_sky >= SKY_EXCELLENT: return 'Excellent'
    if v_sky >= SKY_GOOD:      return 'Good'
    if v_sky >= SKY_FAIR:      return 'Fair'
    return 'Poor'

def lunar_sky_brightness(obs, dt, ra_str, dec_str, dark_sky_mag, k=EXTINCTION_K):
    """
    compute total V-band sky brightness at the target position.
    Krisciunas & Schaefer (1991) PASP 103, 1033.

    returns dict: illum_pct, moon_alt, separation, V_sky, sky_label, moon_up
    """
    obs.date = dt.strftime('%Y/%m/%d %H:%M:%S')

    moon = ephem.Moon()
    moon.compute(obs)
    illum_pct = float(moon.phase)
    moon_alt  = math.degrees(float(moon.alt))

    target = ephem.FixedBody()
    target._ra    = ephem.hours(ra_str)
    target._dec   = ephem.degrees(dec_str)
    target._epoch = ephem.J2000
    target.compute(obs)
    target_alt = math.degrees(float(target.alt))
    sep_deg    = math.degrees(float(ephem.separation(moon, target)))

    # phase angle from illumination fraction
    illum_frac = max(0.0, min(1.0, illum_pct / 100.0))
    alpha_deg  = math.degrees(
        math.acos(max(-1.0, min(1.0, 2.0 * illum_frac - 1.0)))
    )

    V_moon = -12.73 + 0.026 * abs(alpha_deg) + 4e-9 * alpha_deg ** 4
    I_star = 10.0 ** (-0.4 * (V_moon + 16.57))

    rho_deg = max(0.01, sep_deg)
    rho_rad = math.radians(rho_deg)
    f_rho   = (10.0 ** 5.36 * (1.06 + math.cos(rho_rad) ** 2)
               + 10.0 ** (6.15 - rho_deg / 40.0))

    B_dark = 34.08 * 10.0 ** (0.4 * (22.0 - dark_sky_mag))

    X_obj  = _airmass(target_alt)
    X_moon = _airmass(moon_alt)
    moon_up = moon_alt >= 0 and X_moon is not None

    if moon_up and X_obj is not None:
        B_moon = max(0.0, f_rho * I_star
                     * 10.0 ** (-0.4 * k * X_moon)
                     * (1.0 - 10.0 ** (-0.4 * k * X_obj)))
    else:
        B_moon = 0.0

    B_total = B_dark + B_moon
    V_sky   = 22.0 - 2.5 * math.log10(B_total / 34.08)

    return {
        'illum_pct':  illum_pct,
        'moon_alt':   moon_alt,
        'separation': sep_deg,
        'V_sky':      V_sky,
        'sky_label':  _sky_label(V_sky),
        'moon_up':    moon_up,
    }

# ==============================================================================
# CELL 7 — WINDOW CALCULATOR
# ==============================================================================

J2000 = datetime(2000, 1, 1, 12, 0, 0)

def utc_to_jd(dt):
    return 2451545.0 + (dt - J2000).total_seconds() / 86400.0

def gap_occurrences(date):
    """both gap occurrence datetimes for a given UTC date (midnight = start)"""
    gap_dur_h = (GAP_END_PH - GAP_START_PH) * PERIOD_H
    ph_mn     = ((utc_to_jd(date) - T0_JD) * 24.0 / PERIOD_H) % 1.0
    dp        = (GAP_START_PH - ph_mn) % 1.0
    g1s = date + timedelta(hours=dp * PERIOD_H)
    g1e = g1s  + timedelta(hours=gap_dur_h)
    g2s = g1s  + timedelta(hours=PERIOD_H)
    g2e = g1e  + timedelta(hours=PERIOD_H)
    return g1s, g1e, g2s, g2e

def step_through(obs_site, win_s, win_e, radec_lookup):
    if win_s >= win_e:
        return []
    good = []
    t = win_s
    while t <= win_e:
        ra_str, dec_str = interpolate_radec(radec_lookup, t)
        alt = get_alt_ephem(obs_site, t, ra_str, dec_str)
        if alt >= MIN_ALT_DEG:
            good.append((t, alt))
        t += timedelta(minutes=ALT_STEP_MINUTES)
    return good

def compute_window(obs_site, date, cfg, radec_lookup,
                   dark_sky_mag, gap_s=None, gap_e=None):
    """
    compute observable window for one night + site.
    window = dark sky ∩ altitude > MIN_ALT_DEG [∩ phase gap if provided].
    moon metrics computed at window start and midpoint.
    """
    dark_s, dark_e = dark_window(obs_site, date, gap_s)

    win_s = max(dark_s, gap_s) if gap_s is not None else dark_s
    win_e = min(dark_e, gap_e) if gap_e is not None else dark_e

    good = step_through(obs_site, win_s, win_e, radec_lookup)

    result = {
        'date':       date,
        'obs_start':  None,
        'obs_end':    None,
        'minutes':    0,
        'alt_start':  None,
        'alt_end':    None,
        'limit':      '',
        'moon_start': None,
        'moon_mid':   None,
        '_dawn':      dark_e,
        '_gap_e':     gap_e,
    }

    if good:
        t_s, a_s = good[0]
        t_e, a_e = good[-1]
        t_mid    = t_s + (t_e - t_s) / 2
        ra_s,  dec_s = interpolate_radec(radec_lookup, t_s)
        ra_m,  dec_m = interpolate_radec(radec_lookup, t_mid)
        result.update({
            'obs_start':  t_s,
            'obs_end':    t_e,
            'minutes':    len(good) * ALT_STEP_MINUTES,
            'alt_start':  a_s,
            'alt_end':    a_e,
            'limit':      _limit(t_e, gap_e, dark_e),
            'moon_start': lunar_sky_brightness(
                              obs_site, t_s,   ra_s, dec_s, dark_sky_mag),
            'moon_mid':   lunar_sky_brightness(
                              obs_site, t_mid, ra_m, dec_m, dark_sky_mag),
        })
    else:
        if win_s >= win_e:
            result['limit'] = 'not dark'
        else:
            mid = win_s + (win_e - win_s) / 2
            ra_str, dec_str = interpolate_radec(radec_lookup, mid)
            alt_mid = get_alt_ephem(obs_site, mid, ra_str, dec_str)
            result['limit'] = (f"alt {alt_mid:.0f}° at midpoint"
                               if alt_mid < MIN_ALT_DEG else 'no overlap')
    return result

def _limit(obs_end, gap_end, dawn):
    tol = timedelta(minutes=2)
    if gap_end and abs(obs_end - gap_end) < tol: return 'gap ends'
    if abs(obs_end - dawn)                 < tol: return 'dawn'
    return f"alt < {MIN_ALT_DEG:.0f}° at {obs_end.strftime('%H:%M')}"

# ==============================================================================
# CELL 8 — DISPLAY AND CSV
# ==============================================================================

W = 92   # console width

def hm(dt):
    return dt.strftime('%H:%M') if dt else '—'

def hm_local(dt, offset):
    """format time as HH:MM local, noting date change if needed"""
    if dt is None:
        return '—'
    local = utc_to_local(dt, offset)
    return local.strftime('%H:%M')

def rule(char='=', width=W):
    print(char * width)

def section(title, char='='):
    rule(char)
    print(f"  {title}")
    rule(char)

def subsection(title):
    pad = W - 4 - len(title)
    print(f"\n  -- {title} " + '-' * max(pad, 2))

def _moon_line(ms, mm):
    """compact moon sub-line for a result row"""
    if ms is None and mm is None:
        return None
    illum = ms['illum_pct'] if ms else mm['illum_pct']

    def half(m, label):
        if m is None:
            return f"  {label}: —"
        alt_str = f"{m['moon_alt']:>4.0f}° up" if m['moon_up'] else "      down"
        return (f"  {label}: {m['V_sky']:>5.2f} mag/\"  "
                f"{m['sky_label']:<9}"
                f"sep {m['separation']:>4.0f}°  "
                f"moon {alt_str}")

    return (f"           Moon illum {illum:>3.0f}%"
            + half(ms, 'start')
            + half(mm, 'mid'))

# --- ephemeris display -------------------------------------------------------

def print_ephem_table(eph_table, step):
    preferred_cols = ['RA', 'Dec', 'V', 'Delta', 'r', 'Altitude', 'Azimuth']
    cols           = [c for c in preferred_cols if c in eph_table.colnames]
    col_w          = 10

    step_lower = step.lower()
    if   step_lower.endswith('d'):   step_h = float(step_lower[:-1]) * 24
    elif step_lower.endswith('h'):   step_h = float(step_lower[:-1])
    elif step_lower.endswith('min'): step_h = float(step_lower[:-3]) / 60
    else:                            step_h = 24.0

    rows_per_day = max(1, int(round(24.0 / step_h)))

    units = {'RA': '(deg)', 'Dec': '(deg)', 'V': '(mag)',
             'Delta': '(AU)', 'r': '(AU)',
             'Altitude': '(deg)', 'Azimuth': '(deg)'}

    print(f"  {'Date (UTC)':<22}" +
          "".join(f"{c+' '+units.get(c,''):>{col_w}}" for c in cols))
    print("  " + "-" * (22 + col_w * len(cols)))

    for i in range(0, len(eph_table), rows_per_day):
        row  = eph_table[i]
        line = f"  {str(row['Date']):<22}"
        for c in cols:
            try:    line += f"{float(row[c]):>{col_w}.4f}"
            except: line += f"{str(row[c]):>{col_w}}"
        print(line)

    if rows_per_day > 1:
        print(f"\n  note: {len(eph_table)} rows queried ({step} step); "
              f"table shows one row per day. full data → {CSV_EPHEM}")

def save_ephem_csv(eph_table, filename):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([f'# ephemeris: {ASTEROID_ID} {ASTEROID_NAME}',
                         f'step: {EPHEM_STEP}',
                         f'site: {PRIMARY_SITE} ({SITES[PRIMARY_SITE]["mpc_code"]})',
                         f'generated: {datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")}'])
        writer.writerow(eph_table.colnames)
        for row in eph_table:
            writer.writerow([str(row[c]) for c in eph_table.colnames])
    print(f"  saved → {filename}  ({len(eph_table)} rows, "
          f"{len(eph_table.colnames)} columns)")

# --- window display ----------------------------------------------------------

def print_window_row(date, result, cfg, gap_window=''):
    """print window row with UTC and local times, plus moon sub-line"""
    offset = cfg['time_offset']
    mins   = result['minutes']

    if mins > 0:
        utc_rng   = f"{hm(result['obs_start'])} – {hm(result['obs_end'])} UTC"
        local_rng = (f"{hm_local(result['obs_start'], offset)} – "
                     f"{hm_local(result['obs_end'],   offset)} "
                     f"local ({offset_str(offset)})")
        alt_rng   = (f"{result['alt_start']:.0f}° → {result['alt_end']:.0f}°"
                     if result['alt_start'] is not None else '—')
        lim       = result['limit']

        if gap_window:
            print(f"  {date.strftime('%m-%d'):<7} "
                  f"gap {gap_window:<13} "
                  f"{mins:>4} min  "
                  f"{utc_rng:<22}  "
                  f"{local_rng:<26}  "
                  f"{alt_rng:<13}  {lim}")
        else:
            print(f"  {date.strftime('%m-%d'):<7} "
                  f"{mins:>4} min  "
                  f"{utc_rng:<22}  "
                  f"{local_rng:<26}  "
                  f"{alt_rng:<13}  {lim}")

        moon_ln = _moon_line(result.get('moon_start'), result.get('moon_mid'))
        if moon_ln:
            print(moon_ln)
    else:
        lim = result['limit']
        if gap_window:
            print(f"  {date.strftime('%m-%d'):<7} "
                  f"gap {gap_window:<13}    0 min  —  ({lim})")
        else:
            print(f"  {date.strftime('%m-%d'):<7}    0 min  —  ({lim})")

def print_window_header(with_gap=False):
    gap_col = f"{'Gap window':<17}" if with_gap else ''
    print(f"  {'Date':<7} {gap_col}"
          f"{'Min':>4}      "
          f"{'UTC window':<22}  "
          f"{'Local window':<26}  "
          f"{'Altitude':<13}  Limit")
    print("  " + "-" * (W - 2))

# --- CSV for windows ---------------------------------------------------------

def _moon_csv(m, prefix):
    if m is None:
        return {f'{prefix}{k}': '' for k in
                ['v_sky', 'sky_label', 'illum_pct', 'moon_alt', 'separation']}
    return {
        f'{prefix}v_sky':      f"{m['V_sky']:.2f}",
        f'{prefix}sky_label':  m['sky_label'],
        f'{prefix}illum_pct':  f"{m['illum_pct']:.1f}",
        f'{prefix}moon_alt':   f"{m['moon_alt']:.1f}",
        f'{prefix}separation': f"{m['separation']:.1f}",
    }

def _window_row_dict(site, date, result, cfg, gap_label='', gap_window=''):
    offset = cfg['time_offset']
    row = {
        'site':            site,
        'date':            date.strftime('%Y-%m-%d'),
        'gap_label':       gap_label,
        'gap_window':      gap_window,
        'obs_start_utc':   hm(result['obs_start']),
        'obs_end_utc':     hm(result['obs_end']),
        'obs_start_local': hm_local(result['obs_start'], offset),
        'obs_end_local':   hm_local(result['obs_end'],   offset),
        'time_offset':     offset_str(offset),
        'duration_min':    result['minutes'],
        'alt_start':       f"{result['alt_start']:.1f}" if result['alt_start'] is not None else '',
        'alt_end':         f"{result['alt_end']:.1f}"   if result['alt_end']   is not None else '',
        'limit':           result['limit'],
    }
    row.update(_moon_csv(result.get('moon_start'), 'start_'))
    row.update(_moon_csv(result.get('moon_mid'),   'mid_'))
    return row

def save_windows_csv(rows, filename):
    if not rows:
        return
    fields = [
        'site', 'date', 'gap_label', 'gap_window',
        'obs_start_utc', 'obs_end_utc',
        'obs_start_local', 'obs_end_local', 'time_offset',
        'duration_min', 'alt_start', 'alt_end', 'limit',
        'start_v_sky', 'start_sky_label', 'start_illum_pct',
        'start_moon_alt', 'start_separation',
        'mid_v_sky',   'mid_sky_label',   'mid_illum_pct',
        'mid_moon_alt',   'mid_separation',
    ]
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print(f"  saved → {filename}  ({len(rows)} rows)")

# ==============================================================================
# CELL 9 — MAIN
# ==============================================================================

def main():

    # ── run header ────────────────────────────────────────────────────────────
    rule()
    print(f"  ASTEROID OBSERVING SCHEDULER")
    print(f"  Target   : {ASTEROID_ID} {ASTEROID_NAME}")
    print(f"  Dates    : {DATE_START.strftime('%Y-%m-%d')}  to  "
          f"{(DATE_START + timedelta(days=N_DAYS-1)).strftime('%Y-%m-%d')}  "
          f"({N_DAYS} nights)")
    if USE_LIGHTCURVE:
        gap_dur_h     = (GAP_END_PH - GAP_START_PH) * PERIOD_H
        drift_min_day = ((24.0 / PERIOD_H) % 1.0) * PERIOD_H * 60
        print(f"  Mode     : phase gap filter ON")
        print(f"  T0       : JD {T0_JD}   Period: {PERIOD_H} h")
        print(f"  Gap      : phase {GAP_START_PH} – {GAP_END_PH}  "
              f"({gap_dur_h * 60:.1f} min = "
              f"{int(gap_dur_h)}h {int((gap_dur_h % 1) * 60)}m)  "
              f"drift {drift_min_day:.2f} min/day earlier")
    else:
        print(f"  Mode     : full visibility window (no phase filter)")
    print(f"  Alt min  : > {MIN_ALT_DEG:.0f}°   step: {ALT_STEP_MINUTES} min")
    print(f"  Sky model: Krisciunas-Schaefer (1991)   k = {EXTINCTION_K}")
    print(f"  Quality  : Excellent ≥ {SKY_EXCELLENT}  "
          f"Good ≥ {SKY_GOOD}  "
          f"Fair ≥ {SKY_FAIR}  "
          f"Poor < {SKY_FAIR}  (mag/arcsec²)")
    rule()

    # ── ephemeris ─────────────────────────────────────────────────────────────
    print()
    primary_cfg = SITES[PRIMARY_SITE]
    eph_table   = query_ephemeris(
        ASTEROID_ID, primary_cfg['mpc_code'],
        DATE_START, N_DAYS, EPHEM_STEP
    )

    print()
    section(f"EPHEMERIS  ·  {ASTEROID_ID} {ASTEROID_NAME}  ·  "
            f"{PRIMARY_SITE} ({primary_cfg['mpc_code']})  ·  step: {EPHEM_STEP}")
    print()
    print_ephem_table(eph_table, EPHEM_STEP)
    print()
    if SAVE_CSV:
        save_ephem_csv(eph_table, CSV_EPHEM)

    radec_lookup = build_radec_lookup(eph_table)

    # ── observing windows ─────────────────────────────────────────────────────
    all_rows = []

    for site_name, cfg in SITES.items():
        obs = make_obs(cfg)

        site_dark = (DARK_SKY_MAG.get(site_name, 21.5)
                     if isinstance(DARK_SKY_MAG, dict)
                     else float(DARK_SKY_MAG))

        offset     = cfg['time_offset']
        lat_str    = f"{abs(cfg['lat']):.4f}° {'N' if cfg['lat'] >= 0 else 'S'}"
        lon_str    = f"{abs(cfg['lon']):.4f}° {'E' if cfg['lon'] >= 0 else 'W'}"

        print()
        section(f"SITE: {site_name.upper()}")
        print(f"  Location  : {lat_str}   {lon_str}   {cfg['elev']} m")
        print(f"  MPC code  : {cfg['mpc_code']}")
        print(f"  Time zone : {offset_str(offset)}")
        print(f"  Dark sky  : {site_dark} mag/arcsec²")
        rule('-')

        if USE_LIGHTCURVE:
            for gap_label, use_gap2 in [("GAP 1", False), ("GAP 2", True)]:
                subsection(gap_label)
                print()
                print_window_header(with_gap=True)

                for day_offset in range(N_DAYS):
                    d = DATE_START + timedelta(days=day_offset)
                    g1s, g1e, g2s, g2e = gap_occurrences(d)
                    gs, ge = (g2s, g2e) if use_gap2 else (g1s, g1e)

                    result      = compute_window(
                        obs, d, cfg, radec_lookup, site_dark,
                        gap_s=gs, gap_e=ge
                    )
                    gap_win_str = f"{hm(gs)}-{hm(ge)}"
                    print_window_row(d, result, cfg, gap_window=gap_win_str)
                    all_rows.append(
                        _window_row_dict(site_name, d, result, cfg,
                                         gap_label, gap_win_str)
                    )
        else:
            subsection("VISIBILITY WINDOW")
            print()
            print_window_header(with_gap=False)

            for day_offset in range(N_DAYS):
                d      = DATE_START + timedelta(days=day_offset)
                result = compute_window(obs, d, cfg, radec_lookup, site_dark)
                print_window_row(d, result, cfg)
                all_rows.append(
                    _window_row_dict(site_name, d, result, cfg)
                )

    # ── save windows CSV ──────────────────────────────────────────────────────
    print()
    rule()
    if SAVE_CSV and all_rows:
        save_windows_csv(all_rows, CSV_WINDOWS)
    rule()


if __name__ == '__main__':
    main()
