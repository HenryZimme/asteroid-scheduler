# asteroid-scheduler

**Phase-gap-aware multi-site observing window scheduler for asteroid rotational lightcurve campaigns.**

<!-- Published to the Astrophysics Source Code Library: [ascl:XXXX.XXX] -->

## The problem

Single-site lightcurve campaigns on slowly-rotating asteroids frequently leave phase gaps: intervals of the rotation period that are perpetually unobservable from a fixed longitude. This happens because the asteroid's synodic drift rate and the local night window interact to exclude the same phase range night after night. Closing a gap requires coordinating observations from geographically separated sites during the specific UTC windows when the gap phase is actually visible.

No existing tool in the Python small-body astronomy ecosystem addressed this directly. I built this one.

## What it does

Given a known or candidate rotation period, the scheduler identifies when and where each phase gap is observable across any set of ground-based sites. It pulls live ephemerides from the Minor Planet Center, computes per-night dark windows using dynamically calculated astronomical twilight, and evaluates lunar sky brightness using the Krisciunas-Schaefer (1991) physical model. It also runs without lightcurve parameters as a general multi-site visibility calculator.

Output is printed to the console and saved to CSV, with one row per site, night, and gap occurrence. Each row includes the observable window in UTC and local time, asteroid altitude, and full sky brightness metrics (V-band mag/arcsec², moon illumination, separation, and altitude at window start and midpoint).

## Validation: 7605 Cindygraber

Developed for an ongoing photometric campaign on main-belt asteroid 7605 Cindygraber (H = 13.6) from Phillips Academy Observatory (MPC I12). With a rotation period of 11.939 hours and a phase gap spanning 0.72 to 1.00 (about 200 minutes), PAO alone cannot close the gap. The scheduler identified the Canary Islands (G40) and Australia (E62) as ideal complementary sites covering both daily gap occurrences across all 10 nights of the March 2026 campaign window.

The campaign results are being prepared for submission to the Minor Planet Bulletin.

## Visualization

You can find a helpful visualization of the asteroid's observational window from Andover on my website: [asteroid_observer](https://henryzimmerman.net/asteroid_observer).


## Citation

Zimmerman, H. (2026). *asteroid-scheduler: Phase-gap-aware multi-site observing window scheduler for asteroid rotational lightcurve campaigns*. Astrophysics Source Code Library. [ascl:XXXX.XXX]

A BibTeX entry is in `CITATION.cff`.

## Author

Henry Zimmerman  
Phillips Academy, Andover, MA  
MPC Observatory Code: I12
