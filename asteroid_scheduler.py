"""
asteroid ephemeris + observing window scheduler
with Krisciunas-Schaefer (1991) lunar sky brightness model
==========================================================
queries MPC via astroquery.mpc.MPC.get_ephemeris, computes visibility
windows at each site (twilight, altitude), optionally constrains to a
lightcurve phase gap, and reports lunar sky brightness at window start
and midpoint. local times are shown alongside UTC using zoneinfo (stdlib)
with the tzdata package for DST handling — no manual offset updates needed.

phase calculation uses BJD_TDB with light travel time correction (LTTC)
per Eastman et al. (2010). T0_JD must be in BJD_TDB+LTTC format.

install: pip install astroquery astropy ephem tzdata
"""

# ==============================================================================
# CELL 1 — USER CONFIGURATION (only cell that changes between asteroids)
# ==============================================================================

# --- target ---
# ex 1: standard main belt (7605 Cindygraber)
'''ASTEROID_ID   = '7605'
ASTEROID_NAME = 'Cindygraber'
T0_JD        = 2461084.655574
PERIOD_H     = 11.939
PERIOD_ERR_H = 0.001  # period uncertainty in hours #'''

#ex 2: fast rotator / faint target (uncomment to test)
'''ASTEROID_ID   = '1998 KY26'
ASTEROID_NAME = '1998 KY26'
T0_JD        = 2451000.0
PERIOD_H     = 0.0891933333  # 5.3 minutes!
PERIOD_ERR_H = 0.000001#'''

# ex 3: bright / slow rotator (uncomment to test)
ASTEROID_ID   = '343'
ASTEROID_NAME = 'Ostara'

# !! IMPORTANT — BJD_TDB + LTTC EPOCH REQUIRED !!
# T0_JD must be expressed in Barycentric Julian Date (BJD_TDB) with the
# one-way asteroid-to-earth light travel time already subtracted — i.e. the
# epoch represents the moment the asteroid *emitted* the reference light, not
# when it was recorded at your telescope.
#
# if your T0 came from a standard observatory reduction (JD_UTC or JD_TDB
# without LTTC), convert it before entering it here. use the LTTC calculator
# at the bottom of this cell or tools such as:
#   https://astroutils.astronomy.osu.edu/time/bjdconvert.html  (Eastman 2010)
#   or call get_asteroid_emission_bjd() on your original epoch.
#
# failure to use a BJD_TDB+LTTC epoch will introduce a systematic phase
# offset of up to ~8 minutes (delta_au * 0.0057755 days) that drifts with
# changing earth-asteroid distance.
T0_JD        = 2452900.0   # <-- must be BJD_TDB+LTTC; see warning above
PERIOD_H     = 109.9
PERIOD_ERR_H = 0.52

# --- date range ---
from datetime import datetime, timedelta, timezone
# CRITICAL: use UTC midnight, not local time.
# datetime.now() returns naive local time; if passed to astropy Time(..., scale='utc')
# it is silently misinterpreted as UTC, introducing a systematic phase offset equal
# to the local UTC offset (e.g. 4h in EDT → ~0.33 phase error for P=11.939h Cindygraber).
# Normalising to UTC midnight also removes the time-of-day component so that the
# phase reference point is consistent across all sites and all nights.
DATE_START = datetime.now(timezone.utc).replace(
    hour=0, minute=0, second=0, microsecond=0, tzinfo=None
)
N_DAYS     = 15

# --- ephemeris sampling ---
# step passed to MPC.get_ephemeris. examples: '1d', '6h', '1h', '30min'
# '1h' recommended — finer steps give better RA/Dec/Delta interpolation.
EPHEM_STEP = '1h'

# resolution for altitude/moon stepping in the window calculator (minutes)
ALT_STEP_MINUTES = 1

# --- observing sites ---
# tz: IANA timezone string — DST is handled automatically, no seasonal edits needed.
SITES = {
    'G40 Canary 1': {
        'lat':      28.29970,
        'lon':     -16.5082,
        'elev':     2390,
        'mpc_code': 'G40',
    },
    'W88 Chile 2': {
        'lat':      -33.269,
        'lon':      -70.53,
        'elev':     1492,
        'mpc_code': 'W88',
    },
    'E62 Australia 1': {
        'lat':      -31.2816709,
        'lon':      149.0801825,
        'elev':     805,
        'mpc_code': 'E62',
    },
    'I12 Phillips Academy Observatory': {
        'lat':      42.647611,
        'lon':     -71.129,
        'elev':     80,
        'mpc_code': 'I12',
    },
}

# site used for the MPC ephemeris query and for the barycentric correction
# in gap_occurrences. any valid MPC site works; primary site is conventional.
PRIMARY_SITE = 'G40 Canary 1'

# --- altitude threshold ---
MIN_ALT_DEG = 20.0

# --- moon / sky brightness (Krisciunas-Schaefer 1991) ---
DARK_SKY_MAG = {
    'G40 Canary 1':                     21.9,
    'E62 Australia 1':                  21.7,
    'W88 Chile 2':                      21.8,
    'I12 Phillips Academy Observatory': 20.5,
}

EXTINCTION_K = {
    'G40 Canary 1':                     0.15,
    'W88 Chile 2':                      0.15,
    'E62 Australia 1':                  0.20,
    'I12 Phillips Academy Observatory': 0.25,
}

# --- lightcurve parameters ---
USE_LIGHTCURVE = True

GAP_START_PH = 0.72
GAP_END_PH   = 1.00

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
from zoneinfo import ZoneInfo
from astroquery.mpc import MPC
import pytz
from icalendar import Calendar, Event
from timezonefinder import TimezoneFinder
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# astropy: required for BJD_TDB conversion and barycentric light travel time
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

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
    return (dark_start, dark_end) = astronomical twilight end to dawn.
    uses gap_s to determine which night's dark window to use:
      gap in utc morning (hour < 12) -> previous evening's twilight -> this morning's dawn
      gap in utc evening/night        -> this evening's twilight -> next morning's dawn
    when gap_s is None (no lightcurve), defaults to evening.
    """
    if gap_s is not None and gap_s.hour < 12:
        dark_s = evening_twilight_end(obs, date - timedelta(days=1))
        dark_e = morning_dawn(obs, date)
    else:
        dark_s = evening_twilight_end(obs, date)
        dark_e = morning_dawn(obs, date + timedelta(days=1))
    return dark_s, dark_e

def utc_to_local(dt, tz_str):
    """
    convert a naive UTC datetime to a timezone-aware local datetime.
    DST transitions are handled automatically by pytz.
    """
    utc_dt = pytz.utc.localize(dt)
    return utc_dt.astimezone(pytz.timezone(tz_str))

def offset_str(dt, tz_str):
    """
    return the UTC offset string for a given UTC datetime and IANA tz string.
    reflects actual DST state at that moment, e.g. 'UTC-4' in EDT, 'UTC-5' in EST.
    """
    local = utc_to_local(dt, tz_str)
    offset_secs = int(local.utcoffset().total_seconds())
    hours = offset_secs // 3600
    mins  = abs(offset_secs % 3600) // 60
    if mins:
        sign = '+' if hours >= 0 else '-'
        return f"UTC{sign}{abs(hours)}:{mins:02d}"
    sign = '+' if hours >= 0 else ''
    return f"UTC{sign}{hours}"

tf = TimezoneFinder()

def get_site_tz(lat, lon):
    """auto-lookup IANA timezone string from coordinates."""
    return tf.timezone_at(lat=lat, lng=lon)

# ==============================================================================
# CELL 5 — RA/Dec/Delta INTERPOLATION FROM MPC TABLE
# ==============================================================================

def build_radec_delta_lookup(eph_table):
    """
    build sorted list of (datetime, ra_deg, dec_deg, delta_au) from MPC table.
    delta (earth-asteroid distance in AU) is required for LTTC phase correction.
    """
    lookup = []
    for row in eph_table:
        dt    = row['Date'].to_datetime()
        ra    = float(row['RA'])
        dec   = float(row['Dec'])
        delta = float(row['Delta'])
        lookup.append((dt, ra, dec, delta))
    return sorted(lookup, key=lambda x: x[0])

def _interp_raw(lookup, dt):
    """
    shared linear interpolation core. returns (ra_deg, dec_deg, delta_au).
    handles boundary clamping and walks the sorted lookup list.
    """
    if dt <= lookup[0][0]:
        return lookup[0][1], lookup[0][2], lookup[0][3]
    if dt >= lookup[-1][0]:
        return lookup[-1][1], lookup[-1][2], lookup[-1][3]
    for i in range(len(lookup) - 1):
        t0, r0, d0, dl0 = lookup[i]
        t1, r1, d1, dl1 = lookup[i + 1]
        if t0 <= dt <= t1:
            frac    = (dt - t0).total_seconds() / (t1 - t0).total_seconds()
            ra_deg  = r0  + frac * (r1  - r0)
            dec_deg = d0  + frac * (d1  - d0)
            delta   = dl0 + frac * (dl1 - dl0)
            return ra_deg, dec_deg, delta
    return lookup[-1][1], lookup[-1][2], lookup[-1][3]

def _deg_to_ephem_str(ra_deg, dec_deg):
    """convert decimal ra/dec degrees to ephem-compatible h:m:s / ±d:m:s strings."""
    ra_h   = ra_deg / 15.0
    ra_str = (f"{int(ra_h)}:{int((ra_h % 1) * 60)}:"
              f"{((ra_h % 1) * 60 % 1) * 60:.2f}")
    sign    = '+' if dec_deg >= 0 else '-'
    dec_abs = abs(dec_deg)
    dec_str = (f"{sign}{int(dec_abs)}:{int((dec_abs % 1) * 60)}:"
               f"{((dec_abs % 1) * 60 % 1) * 60:.2f}")
    return ra_str, dec_str

def interpolate_radec(lookup, dt):
    """
    linear interpolation of RA/Dec for a given datetime.
    returns (ra_str, dec_str) in ephem-compatible format.
    used by get_alt_ephem and lunar_sky_brightness (ephem callers).
    """
    ra_deg, dec_deg, _ = _interp_raw(lookup, dt)
    return _deg_to_ephem_str(ra_deg, dec_deg)

def interpolate_position(lookup, dt):
    """
    linear interpolation returning physical coordinates for BJD calculation.
    returns (ra_deg, dec_deg, delta_au).
    """
    return _interp_raw(lookup, dt)

# ==============================================================================
# CELL 5B — BJD_TDB + LTTC PHASE CORRECTION
# ==============================================================================

def get_asteroid_emission_bjd(utc_dt, ra_deg, dec_deg, delta_au, site_cfg=None):
    """
    convert a UTC observation time to BJD_TDB, then subtract the one-way
    asteroid-to-earth light travel time (LTTC) to obtain the epoch at which
    the asteroid *emitted* the observed light.

    this is the correct time standard for comparing against a T0_JD that was
    derived from a phased lightcurve using BJD_TDB+LTTC (Eastman et al. 2010).

    args:
        utc_dt   : naive python datetime in UTC
        ra_deg   : asteroid right ascension in decimal degrees (J2000)
        dec_deg  : asteroid declination in decimal degrees (J2000)
        delta_au : earth-asteroid distance in AU at utc_dt
        site_cfg : dict with keys lat, lon, elev for the observatory.
                   if None, geocenter is used (difference is < 0.02 s).

    returns:
        float: BJD_TDB of asteroid emission event
    """
    t  = Time(utc_dt, format='datetime', scale='utc')
    sc = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)

    if site_cfg is not None:
        loc = EarthLocation(lat=site_cfg['lat']  * u.deg,
                            lon=site_cfg['lon']  * u.deg,
                            height=site_cfg['elev'] * u.m)
    else:
        # geocenter: adequate when no site_cfg is provided
        loc = EarthLocation(lat=0 * u.deg, lon=0 * u.deg, height=0 * u.m)

    # barycentric light travel time: accounts for earth's orbital motion
    # relative to the solar system barycenter (roemer delay + relativistic terms)
    ltt_bary  = t.light_travel_time(sc, 'barycentric', location=loc)
    bjd_tdb   = (t.tdb + ltt_bary).jd

    # subtract one-way light travel time from asteroid to earth:
    # 0.0057755 days/AU = 1 AU / c
    lttc_days = delta_au * 0.0057755
    return bjd_tdb - lttc_days

# ==============================================================================
# CELL 6 — KRISCIUNAS-SCHAEFER (1991) LUNAR SKY BRIGHTNESS MODEL
# ==============================================================================

def _airmass(alt_deg):
    """rozenberg (1966) airmass. returns None if more than 1 deg below horizon."""
    if alt_deg < -1.0:
        return None
    z_rad = math.radians(max(0.0, 90.0 - alt_deg))
    return 1.0 / (math.cos(z_rad) + 0.025 * math.exp(-11.0 * math.cos(z_rad)))


def lunar_sky_brightness(obs, dt, ra_str, dec_str, dark_sky_mag, k=0.172):
    """
    compute total V-band sky brightness at the target position.
    Krisciunas & Schaefer (1991) PASP 103, 1033.

    returns dict: illum_pct, moon_alt, separation, V_sky, moon_up
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
        'moon_up':    moon_up,
    }

# ==============================================================================
# CELL 7 — WINDOW CALCULATOR
# ==============================================================================

def gap_occurrences(date, radec_delta_lookup, site_cfg=None):
    """
    compute the UTC start/end times for gap 1 and gap 2 on a given date,
    using BJD_TDB + LTTC for the phase calculation.

    args:
        date               : python datetime (UTC midnight of the observing date)
        radec_delta_lookup : output of build_radec_delta_lookup()
        site_cfg           : site dict (lat, lon, elev) for barycentric correction.
                             uses primary site in main(); geocenter if None.

    returns:
        g1s_min, g1e_max, g2s_min, g2e_max — UTC datetimes, expanded by
        accumulated period uncertainty since T0.
    """
    # Normalise to UTC midnight — eliminates any time-of-day component that
    # would otherwise produce a systematic phase offset equal to hours-since-midnight.
    date = date.replace(hour=0, minute=0, second=0, microsecond=0)

    # interpolate asteroid position at start of this date for phase reference
    ra_deg, dec_deg, delta_au = interpolate_position(radec_delta_lookup, date)

    # convert to BJD_TDB and subtract LTTC to get asteroid emission epoch
    emission_bjd = get_asteroid_emission_bjd(
        date, ra_deg, dec_deg, delta_au, site_cfg
    )

    # phase and accumulated period drift since T0
    days_since_t0       = emission_bjd - T0_JD
    cycles              = (days_since_t0 * 24.0) / PERIOD_H
    accumulated_error_h = abs(cycles) * PERIOD_ERR_H

    gap_dur_h = (GAP_END_PH - GAP_START_PH) * PERIOD_H
    ph_mn     = (days_since_t0 * 24.0 / PERIOD_H) % 1.0
    dp        = (GAP_START_PH - ph_mn) % 1.0

    # nominal gap windows (UTC)
    g1s = date + timedelta(hours=dp * PERIOD_H)
    g1e = g1s + timedelta(hours=gap_dur_h)
    g2s = g1s + timedelta(hours=PERIOD_H)
    g2e = g1e + timedelta(hours=PERIOD_H)

    # expand by accumulated drift uncertainty so observer doesn't miss the gap
    g1s_min = g1s - timedelta(hours=accumulated_error_h)
    g1e_max = g1e + timedelta(hours=accumulated_error_h)
    g2s_min = g2s - timedelta(hours=accumulated_error_h)
    g2e_max = g2e + timedelta(hours=accumulated_error_h)

    return g1s_min, g1e_max, g2s_min, g2e_max

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
        k_val = EXTINCTION_K.get(cfg.get('name', ''), 0.172)
        result.update({
            'obs_start':  t_s,
            'obs_end':    t_e,
            'minutes':    len(good) * ALT_STEP_MINUTES,
            'alt_start':  a_s,
            'alt_end':    a_e,
            'limit':      _limit(t_e, gap_e, dark_e),
            'moon_start': lunar_sky_brightness(
                              obs_site, t_s,   ra_s, dec_s, dark_sky_mag, k=k_val),
            'moon_mid':   lunar_sky_brightness(
                              obs_site, t_mid, ra_m, dec_m, dark_sky_mag, k=k_val),
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

def hm_local(dt, tz_str):
    """
    format a naive UTC datetime as HH:MM in the site's local timezone.
    DST offset is computed at the actual moment dt, not a fixed integer.
    """
    if dt is None:
        return '—'
    return utc_to_local(dt, tz_str).strftime('%H:%M')

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
              f"table shows one row per day. full data -> {CSV_EPHEM}")

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
    print(f"  saved -> {filename}  ({len(eph_table)} rows, "
          f"{len(eph_table.colnames)} columns)")

# --- window display ----------------------------------------------------------

def print_window_row(date, result, cfg, gap_window=''):
    """print window row with UTC and local times, plus moon sub-line"""
    tz_str = cfg['tz']
    mins   = result['minutes']

    if mins > 0:
        tz_label  = offset_str(result['obs_start'], tz_str)
        utc_rng   = f"{hm(result['obs_start'])} – {hm(result['obs_end'])} UTC"
        local_rng = (f"{hm_local(result['obs_start'], tz_str)} – "
                     f"{hm_local(result['obs_end'],   tz_str)} "
                     f"local ({tz_label})")
        alt_rng   = (f"{result['alt_start']:.0f}° -> {result['alt_end']:.0f}°"
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
                ['v_sky', 'illum_pct', 'moon_alt', 'separation']}
    return {
        f'{prefix}v_sky':      f"{m['V_sky']:.2f}",
        f'{prefix}illum_pct':  f"{m['illum_pct']:.1f}",
        f'{prefix}moon_alt':   f"{m['moon_alt']:.1f}",
        f'{prefix}separation': f"{m['separation']:.1f}",
    }

def _window_row_dict(site, date, result, cfg, gap_label='', gap_window=''):
    tz_str = cfg['tz']
    ref_dt = result['obs_start'] if result['obs_start'] else date
    tz_label = offset_str(ref_dt, tz_str)
    row = {
        'site':            site,
        'date':            date.strftime('%Y-%m-%d'),
        'gap_label':       gap_label,
        'gap_window':      gap_window,
        'obs_start_utc':   hm(result['obs_start']),
        'obs_end_utc':     hm(result['obs_end']),
        'obs_start_local': hm_local(result['obs_start'], tz_str),
        'obs_end_local':   hm_local(result['obs_end'],   tz_str),
        'time_offset':     tz_label,
        'duration_min':    result['minutes'],
        'alt_start':       f"{result['alt_start']:.1f}" if result['alt_start'] is not None else '',
        'alt_end':         f"{result['alt_end']:.1f}"   if result['alt_end']   is not None else '',
        'limit':           result['limit'],
    }
    row.update(_moon_csv(result.get('moon_start'), 'start_'))
    row.update(_moon_csv(result.get('moon_mid'), 'mid_'))

    # preserve raw datetime objects for plotting and calendar export
    row['_obs_start_dt'] = result['obs_start']
    row['_obs_end_dt']   = result['obs_end']

    return row

def save_windows_csv(rows, filename):
    if not rows:
        return
    fields = [
        'site', 'date', 'gap_label', 'gap_window',
        'obs_start_utc', 'obs_end_utc',
        'obs_start_local', 'obs_end_local', 'time_offset',
        'duration_min', 'alt_start', 'alt_end', 'limit',
        'start_v_sky',  'start_illum_pct',
        'start_moon_alt', 'start_separation',
        'mid_v_sky', 'mid_illum_pct',
        'mid_moon_alt',   'mid_separation',
    ]
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print(f"  saved -> {filename}  ({len(rows)} rows)")

def save_windows_ics(rows, filename_prefix):
    """exports observable windows to .ics calendar files per site."""
    if not rows: return
    cal_per_site = {}

    for row in rows:
        site = row['site']
        if site not in cal_per_site:
            cal_per_site[site] = Calendar()
            cal_per_site[site].add('prodid', '-//Asteroid Observing Scheduler//')
            cal_per_site[site].add('version', '2.0')

        start_dt = row.get('_obs_start_dt')
        end_dt   = row.get('_obs_end_dt')

        if start_dt and end_dt:
            e = Event()
            e.add('summary', f"Observe {ASTEROID_NAME} - {row['gap_label']}")
            e.add('dtstart', pytz.utc.localize(start_dt))
            e.add('dtend',   pytz.utc.localize(end_dt))
            description = (f"Moon Illumination: {row.get('start_illum_pct')}%\n"
                           f"Moon Altitude: {row.get('start_moon_alt')}°\n"
                           f"Sky Brightness: {row.get('start_v_sky')} mag/arcsec²\n"
                           f"Limit Constraint: {row['limit']}")
            e.add('description', description)
            cal_per_site[site].add_component(e)

    for site, cal in cal_per_site.items():
        safe_site = site.replace(' ', '_').replace('/', '_')
        fname = f"{filename_prefix}_{safe_site}.ics"
        with open(fname, 'wb') as f:
            f.write(cal.to_ical())
        print(f"  saved calendar: {fname}")

def plot_observing_timeline(rows):
    """gantt-style timeline chart of observing windows."""
    if not rows: return
    fig, ax = plt.subplots(figsize=(10, 4))
    sites  = list(set(r['site'] for r in rows))
    site_y = {site: i for i, site in enumerate(sites)}

    for row in rows:
        start_dt = row.get('_obs_start_dt')
        end_dt   = row.get('_obs_end_dt')
        if start_dt and end_dt:
            duration_days = (end_dt - start_dt).total_seconds() / 86400
            ax.barh(site_y[row['site']], duration_days,
                    left=start_dt, height=0.4, color='royalblue')

    ax.set_yticks(range(len(sites)))
    ax.set_yticklabels(sites)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
    plt.xticks(rotation=45)
    plt.xlabel('UTC Time')
    plt.title(f'Observing Windows Timeline: {ASTEROID_NAME}')
    plt.grid(axis='x', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(f'timeline_{ASTEROID_ID}.png')
    plt.show()

def plot_sky_brightness(rows):
    """plots lunar sky brightness at the start of each window."""
    if not rows: return
    fig, ax = plt.subplots(figsize=(10, 4))
    sites = list(set(r['site'] for r in rows))

    for site in sites:
        site_rows = [r for r in rows
                     if r['site'] == site and r.get('_obs_start_dt') is not None]
        if not site_rows: continue
        times = [r['_obs_start_dt'] for r in site_rows]
        v_sky = [float(r['start_v_sky']) for r in site_rows if r['start_v_sky']]
        ax.plot(times, v_sky, marker='o', linestyle='-', label=site)

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    plt.xticks(rotation=45)
    plt.ylabel('Sky Brightness (V mag/arcsec²)')
    plt.title(f'Sky Brightness at Window Start: {ASTEROID_NAME}')
    plt.legend()
    plt.grid(linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(f'sky_brightness_{ASTEROID_ID}.png')
    plt.show()

# ==============================================================================
# CELL 9 — MAIN
# ==============================================================================

def main():

    # -- run header ------------------------------------------------------------
    rule()
    print(f"  ASTEROID OBSERVING SCHEDULER")
    print(f"  Target   : {ASTEROID_ID} {ASTEROID_NAME}")
    print(f"  Dates    : {DATE_START.strftime('%Y-%m-%d')}  to  "
          f"{(DATE_START + timedelta(days=N_DAYS-1)).strftime('%Y-%m-%d')}  "
          f"({N_DAYS} nights)")
    print(f"  Phase    : BJD_TDB + LTTC  (T0 must be in same standard)")
    if USE_LIGHTCURVE:
        gap_dur_h     = (GAP_END_PH - GAP_START_PH) * PERIOD_H
        drift_min_day = ((24.0 / PERIOD_H) % 1.0) * PERIOD_H * 60
        print(f"  Mode     : phase gap filter ON")
        print(f"  T0       : BJD_TDB {T0_JD}   Period: {PERIOD_H} h")
        print(f"  Gap      : phase {GAP_START_PH} – {GAP_END_PH}  "
              f"({gap_dur_h * 60:.1f} min = "
              f"{int(gap_dur_h)}h {int((gap_dur_h % 1) * 60)}m)  "
              f"drift {drift_min_day:.2f} min/day earlier")
    else:
        print(f"  Mode     : full visibility window (no phase filter)")
    print(f"  Alt min  : > {MIN_ALT_DEG:.0f}°   step: {ALT_STEP_MINUTES} min")
    print(f"  Sky model: Krisciunas-Schaefer (1991)   k = {EXTINCTION_K}")
    rule()

    # -- ephemeris -------------------------------------------------------------
    print()
    primary_cfg = SITES[PRIMARY_SITE]
    eph_table   = query_ephemeris(
        ASTEROID_ID, primary_cfg['mpc_code'],
        DATE_START, N_DAYS, EPHEM_STEP
    )

    print()
    section(f"EPHEMERIS  .  {ASTEROID_ID} {ASTEROID_NAME}  .  "
            f"{PRIMARY_SITE} ({primary_cfg['mpc_code']})  .  step: {EPHEM_STEP}")
    print()
    print_ephem_table(eph_table, EPHEM_STEP)
    print()
    if SAVE_CSV:
        save_ephem_csv(eph_table, CSV_EPHEM)

    # build unified lookup: (dt, ra_deg, dec_deg, delta_au)
    # used by both interpolate_radec (ephem strings) and interpolate_position (degrees + delta)
    radec_delta_lookup = build_radec_delta_lookup(eph_table)

    # -- observing windows -----------------------------------------------------
    all_rows = []

    tf_finder = TimezoneFinder()
    for site_name, cfg in SITES.items():
        obs = make_obs(cfg)

        site_dark = (DARK_SKY_MAG.get(site_name, 21.5)
                     if isinstance(DARK_SKY_MAG, dict)
                     else float(DARK_SKY_MAG))

        # auto-lookup timezone if not set
        if 'tz' not in cfg:
            cfg['tz'] = tf_finder.timezone_at(lng=cfg['lon'], lat=cfg['lat'])
        tz_str = cfg['tz']

        lat_str = f"{abs(cfg['lat']):.4f}° {'N' if cfg['lat'] >= 0 else 'S'}"
        lon_str = f"{abs(cfg['lon']):.4f}° {'E' if cfg['lon'] >= 0 else 'W'}"
        tz_label_start = offset_str(DATE_START, tz_str)

        print()
        section(f"SITE: {site_name.upper()}")
        print(f"  Location  : {lat_str}   {lon_str}   {cfg['elev']} m")
        print(f"  MPC code  : {cfg['mpc_code']}")
        print(f"  Time zone : {tz_str}  ({tz_label_start} at campaign start)")
        print(f"  Dark sky  : {site_dark} mag/arcsec²")
        rule('-')

        if USE_LIGHTCURVE:
            for gap_label, use_gap2 in [("GAP 1", False), ("GAP 2", True)]:
                subsection(gap_label)
                print()
                print_window_header(with_gap=True)

                for day_offset in range(N_DAYS):
                    d = DATE_START + timedelta(days=day_offset)

                    # phase uses BJD_TDB+LTTC via the primary site's location
                    g1s, g1e, g2s, g2e = gap_occurrences(
                        d, radec_delta_lookup, primary_cfg
                    )
                    gs, ge = (g2s, g2e) if use_gap2 else (g1s, g1e)

                    result      = compute_window(
                        obs, d, cfg, radec_delta_lookup, site_dark,
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
                result = compute_window(obs, d, cfg, radec_delta_lookup, site_dark)
                print_window_row(d, result, cfg)
                all_rows.append(
                    _window_row_dict(site_name, d, result, cfg)
                )

    # -- save windows CSV ------------------------------------------------------
    print()
    rule()

    if SAVE_CSV:
        save_windows_csv(all_rows, CSV_WINDOWS)
        print(f"\n  >> saved windows to {CSV_WINDOWS}")

    print("\nGenerating actionable outputs...")
    save_windows_ics(all_rows, f"windows_{ASTEROID_ID}")
    plot_observing_timeline(all_rows)
    plot_sky_brightness(all_rows)

    rule()

if __name__ == '__main__':
    main()
