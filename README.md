# lissattracker

A satellite tracker that runs in your terminal with **zero dependencies** — no
skyfield, no sgp4, no numpy, no requests. Just the Python standard library.

SGP4 is implemented from scratch, so the only thing you need is a Python 3.

```sh
pip install git+https://github.com/TOMSPY11/python-sat-tracker
lissattracker --lat 13.82 --lon 100.51 --alt 10
```

> Not on PyPI yet. Once it is published, `pip install lissattracker` will work.

```
 lissattracker  2026-10-23 13:25:57Z  +51.478,-0.002 @45m  2/17 up  [celestrak (stations)]

   SATELLITE            ELEV    AZIM         RANGE   RATE     LAT      LON     ALT
>1 ISS (ZARYA)          42.6°  118.3° ESE     812km  -5.21   46.68    14.02   421km │
 2 CSS (TIANHE)          7.1°  289.4° WNW    1904km  +6.44   51.59   -21.80   387km │
 3 PROGRESS-MS 34       -9.9°  121.9° ESE    3378km  +0.88   20.87    29.38   307km │        N
 4 KNACKSAT-2          -29.4°  162.0° SSE    6959km  -3.13  -11.11    16.26   354km │    .........
                                                                                    │ W     +     E
                                                                                    │        S
 ISS (ZARYA) [25544]  period 92.8min  incl 51.6°  TLE age 0.4d  v 7.66km/s
 j/k move  p passes  m radar  r refetch  space pause  q quit
```

## Usage

```sh
# live TUI (your location is remembered after the first run)
lissattracker --lat 13.82 --lon 100.51 --alt 10

# one-shot table of everything currently above the horizon
lissattracker --once
lissattracker --once --all      # include everything below, too

# when is it next overhead?
lissattracker --passes ISS --hours 48 --min-el 10

# other Celestrak groups, or your own elements
lissattracker --group visual     # also: active, starlink, noaa, gps-ops, ...
lissattracker --tle my_sats.tle
```

Without installing, the module runs straight from a checkout:

```sh
python3 src/lissattracker.py --once
```

### Keys

| key | action |
| --- | --- |
| `j` / `k`, arrows | move the cursor |
| `p` | find the next pass for the selected satellite |
| `m` | toggle the sky radar |
| `r` | refetch TLEs from Celestrak |
| `space` | freeze the clock |
| `q` | quit |

## What it actually does

1. Pulls TLEs from Celestrak (cached on disk for 2 hours, falls back to a stale
   cache if the network is down).
2. Propagates them with SGP4 — un-Kozai'd mean motion, secular J2/J4 terms,
   atmospheric drag, long- and short-period periodics, and a Kepler solve.
3. Rotates TEME into ECEF through GMST, then to WGS84 latitude/longitude/altitude.
4. Converts to topocentric azimuth / elevation / range for your location, plus
   range rate so you can see what is coming toward you.
5. Predicts passes by stepping the horizon coarsely and bisecting each crossing
   down to about a second.

## Accuracy

The propagator is checked against the official SGP4 verification vectors from
Vallado's *Revisiting Spacetrack Report #3* (satellite 00005, five epochs out to
1440 minutes):

| check | result |
| --- | --- |
| position vs. reference | worst error ~7 nanometres |
| velocity vs. reference | worst error < 0.001 mm/s |

Those vectors only validate the TEME state vector, so the reference-frame chain
is checked separately against an independent ISS position service:

| quantity | delta |
| --- | --- |
| ground track | 1.15 km (different TLE epoch) |
| altitude | 110 m |
| speed | < 1 m/s |

## Limits

Deep-space orbits — anything with a period of 225 minutes or more, so GEO,
Molniya, GPS — are **not** propagated. SDP4's lunar-solar resonance terms are
not implemented. Those satellites are reported as unsupported rather than
silently given a wrong position.

Polar motion is ignored in the TEME→ECEF rotation, which is worth a few metres.
Atmospheric refraction near the horizon is not modelled.

## Development

```sh
git clone https://github.com/TOMSPY11/python-sat-tracker
cd python-sat-tracker
PYTHONPATH=src python3 -m unittest discover -s tests -v    # 21 tests, no deps
python3 -m build                                           # sdist + wheel
```

## Licence

MIT
