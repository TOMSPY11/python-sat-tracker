#!/usr/bin/env python3
"""i_dont_wanna_die - a satellite tracker with no dependencies at all.

Pulls TLEs from Celestrak, propagates them with a from-scratch SGP4
(near-Earth model, Spacetrack Report #3 / Vallado), and draws what is over
your head right now as a curses sky-radar plus a live table.

stdlib only: no skyfield, no sgp4, no numpy, no requests.

    python3 i_dont_wanna_die.py --lat 55.75 --lon 37.62 --alt 150
    python3 i_dont_wanna_die.py --group visual
    python3 i_dont_wanna_die.py --once          # one-shot table, no curses
    python3 i_dont_wanna_die.py --passes ISS    # next 48h of passes
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import locale
import math
import os
import sys
import time
import urllib.error
import urllib.request

TAU = 2.0 * math.pi
DEG = 180.0 / math.pi
RAD = math.pi / 180.0

# ---------------------------------------------------------------- WGS72 (SGP4)
RE = 6378.135              # km, earth radius
MU = 398600.8              # km^3/s^2
XKE = 60.0 / math.sqrt(RE * RE * RE / MU)
J2 = 0.001082616
J3 = -0.00000253881
J4 = -0.00000165597
J3OJ2 = J3 / J2
X2O3 = 2.0 / 3.0
VKMPERSEC = RE * XKE / 60.0

# WGS84, used only for turning ECEF into lat/lon/alt
WGS84_A = 6378.137
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)

CACHE_DIR = os.path.expanduser("~/.cache/i_dont_wanna_die")
CONFIG = os.path.expanduser("~/.config/i_dont_wanna_die.json")
CACHE_TTL = 2 * 3600
GP_URL = "https://celestrak.org/NORAD/elements/gp.php?GROUP={}&FORMAT=tle"


# ===========================================================================
# TLE parsing
# ===========================================================================

def _expfield(s: str) -> float:
    """Decode a TLE exponential field like ' 11606-4' -> 0.11606e-4."""
    s = s.strip()
    if not s or s in ("0", "+0", "-0", "00000+0", "00000-0"):
        return 0.0
    sign = -1.0 if s[0] == "-" else 1.0
    s = s.lstrip("+-")
    if "-" in s[1:]:
        mant, _, exp = s.rpartition("-")
        exp = -int(exp)
    elif "+" in s[1:]:
        mant, _, exp = s.rpartition("+")
        exp = int(exp)
    else:
        mant, exp = s, 0
    if not mant:
        return 0.0
    return sign * float("0." + mant) * (10.0 ** exp)


class TLE:
    __slots__ = ("name", "norad", "epoch", "epoch_jd", "bstar", "inclo",
                 "nodeo", "ecco", "argpo", "mo", "no_kozai", "period")

    def __init__(self, name: str, l1: str, l2: str):
        self.name = name.strip() or l2[2:7].strip()
        self.norad = l1[2:7].strip()

        yy = int(l1[18:20])
        doy = float(l1[20:32])
        year = 2000 + yy if yy < 57 else 1900 + yy
        self.epoch = (dt.datetime(year, 1, 1, tzinfo=dt.timezone.utc)
                      + dt.timedelta(days=doy - 1.0))
        self.epoch_jd = julian_day(self.epoch)

        self.bstar = _expfield(l1[53:61])
        self.inclo = float(l2[8:16]) * RAD
        self.nodeo = float(l2[17:25]) * RAD
        self.ecco = float("0." + l2[26:33].strip())
        self.argpo = float(l2[34:42]) * RAD
        self.mo = float(l2[43:51]) * RAD
        rev_per_day = float(l2[52:63])
        self.no_kozai = rev_per_day * TAU / 1440.0   # rad/min
        self.period = 1440.0 / rev_per_day if rev_per_day else 0.0

    @property
    def deep_space(self) -> bool:
        return self.period >= 225.0

    def age_days(self, now: dt.datetime) -> float:
        return (now - self.epoch).total_seconds() / 86400.0


def parse_tle_text(text: str):
    """Accept 3-line (name + 2 element lines) or bare 2-line TLE sets."""
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    out, i = [], 0
    while i < len(lines):
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            name, l1, l2 = "", lines[i], lines[i + 1]
            i += 2
        elif (i + 2 < len(lines) and lines[i + 1].startswith("1 ")
              and lines[i + 2].startswith("2 ")):
            name, l1, l2 = lines[i], lines[i + 1], lines[i + 2]
            i += 3
        else:
            i += 1
            continue
        try:
            out.append(TLE(name, l1, l2))
        except (ValueError, IndexError):
            pass
    return out


# ===========================================================================
# SGP4 - near-earth model
# ===========================================================================

class PropagationError(Exception):
    pass


class SGP4:
    """Near-earth SGP4.  Deep-space (period >= 225 min) is not modelled."""

    def __init__(self, tle: TLE):
        self.tle = tle
        self.deep_space = tle.deep_space
        if self.deep_space:
            return

        e = tle.ecco
        self.inclo = tle.inclo
        self.bstar = tle.bstar

        cosio = math.cos(tle.inclo)
        self.cosio = cosio
        self.sinio = sinio = math.sin(tle.inclo)
        theta2 = cosio * cosio
        theta4 = theta2 * theta2
        self.x3thm1 = x3thm1 = 3.0 * theta2 - 1.0
        self.x1mth2 = x1mth2 = 1.0 - theta2
        self.x7thm1 = 7.0 * theta2 - 1.0

        eosq = e * e
        omeosq = 1.0 - eosq
        rteosq = math.sqrt(omeosq)

        # un-Kozai the mean motion
        no = tle.no_kozai
        a1 = (XKE / no) ** X2O3
        d1 = 0.75 * J2 * x3thm1 / (rteosq * omeosq)
        del_ = d1 / (a1 * a1)
        adel = a1 * (1.0 - del_ * del_ - del_ * (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0))
        del_ = d1 / (adel * adel)
        self.no = no = no / (1.0 + del_)

        ao = (XKE / no) ** X2O3
        self.ao = ao
        po = ao * omeosq
        pinvsq = 1.0 / (po * po)
        rp = ao * (1.0 - e)
        self.perigee_km = (rp - 1.0) * RE

        # drag-term atmosphere fit
        sfour = 78.0 / RE + 1.0
        qzms24 = ((120.0 - 78.0) / RE) ** 4
        if self.perigee_km < 156.0:
            sfour = self.perigee_km - 78.0
            if self.perigee_km < 98.0:
                sfour = 20.0
            qzms24 = ((120.0 - sfour) / RE) ** 4
            sfour = sfour / RE + 1.0

        tsi = 1.0 / (ao - sfour)
        self.eta = eta = ao * e * tsi
        etasq = eta * eta
        eeta = e * eta
        psisq = abs(1.0 - etasq)
        coef = qzms24 * (tsi ** 4)
        coef1 = coef / (psisq ** 3.5)

        cc2 = coef1 * no * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
                            + 0.375 * J2 * tsi / psisq * x3thm1
                            * (8.0 + 3.0 * etasq * (8.0 + etasq)))
        self.cc1 = cc1 = self.bstar * cc2
        cc3 = 0.0
        if e > 1.0e-4:
            cc3 = -2.0 * coef * tsi * J3OJ2 * no * sinio / e
        self.cc4 = 2.0 * no * coef1 * ao * omeosq * (
            eta * (2.0 + 0.5 * etasq) + e * (0.5 + 2.0 * etasq)
            - J2 * tsi / (ao * psisq) * (
                -3.0 * x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                + 0.75 * x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq))
                * math.cos(2.0 * tle.argpo)))
        self.cc5 = 2.0 * coef1 * ao * omeosq * (
            1.0 + 2.75 * (etasq + eeta) + eeta * etasq)

        temp1 = 1.5 * J2 * pinvsq * no
        temp2 = 0.5 * temp1 * J2 * pinvsq
        temp3 = -0.46875 * J4 * pinvsq * pinvsq * no
        self.mdot = (no + 0.5 * temp1 * rteosq * x3thm1
                     + 0.0625 * temp2 * rteosq
                     * (13.0 - 78.0 * theta2 + 137.0 * theta4))
        self.argpdot = (-0.5 * temp1 * (1.0 - 5.0 * theta2)
                        + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4)
                        + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4))
        xhdot1 = -temp1 * cosio
        self.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2)
                                 + 2.0 * temp3 * (3.0 - 7.0 * theta2)) * cosio

        self.omgcof = self.bstar * cc3 * math.cos(tle.argpo)
        self.xmcof = 0.0
        if e > 1.0e-4:
            self.xmcof = -X2O3 * coef * self.bstar / eeta
        self.nodecf = 3.5 * omeosq * xhdot1 * cc1
        self.t2cof = 1.5 * cc1
        if abs(cosio + 1.0) > 1.5e-12:
            self.xlcof = -0.25 * J3OJ2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio)
        else:
            self.xlcof = -0.25 * J3OJ2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12
        self.aycof = -0.5 * J3OJ2 * sinio
        self.delmo = (1.0 + eta * math.cos(tle.mo)) ** 3
        self.sinmao = math.sin(tle.mo)

        # the "simplified" branch, for very low perigees
        self.isimp = self.perigee_km < 220.0
        self.d2 = self.d3 = self.d4 = 0.0
        self.t3cof = self.t4cof = self.t5cof = 0.0
        if not self.isimp:
            cc1sq = cc1 * cc1
            self.d2 = d2 = 4.0 * ao * tsi * cc1sq
            temp = d2 * tsi * cc1 / 3.0
            self.d3 = d3 = (17.0 * ao + sfour) * temp
            self.d4 = d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * cc1
            self.t3cof = d2 + 2.0 * cc1sq
            self.t4cof = 0.25 * (3.0 * d3 + cc1 * (12.0 * d2 + 10.0 * cc1sq))
            self.t5cof = 0.2 * (3.0 * d4 + 12.0 * cc1 * d3 + 6.0 * d2 * d2
                                + 15.0 * cc1sq * (2.0 * d2 + cc1sq))

    # ------------------------------------------------------------------
    def propagate(self, tsince: float):
        """tsince in minutes from TLE epoch -> (r, v) in TEME, km and km/s."""
        if self.deep_space:
            raise PropagationError("deep-space orbit (period >= 225 min)")

        tle = self.tle
        t = tsince
        t2 = t * t

        xmdf = tle.mo + self.mdot * t
        argpdf = tle.argpo + self.argpdot * t
        nodedf = tle.nodeo + self.nodedot * t
        argpm, mm = argpdf, xmdf
        nodem = nodedf + self.nodecf * t2

        tempa = 1.0 - self.cc1 * t
        tempe = self.bstar * self.cc4 * t
        templ = self.t2cof * t2

        if not self.isimp:
            delomg = self.omgcof * t
            delm = self.xmcof * (((1.0 + self.eta * math.cos(xmdf)) ** 3) - self.delmo)
            temp = delomg + delm
            mm = xmdf + temp
            argpm = argpdf - temp
            t3 = t2 * t
            t4 = t3 * t
            tempa -= self.d2 * t2 + self.d3 * t3 + self.d4 * t4
            tempe += self.bstar * self.cc5 * (math.sin(mm) - self.sinmao)
            templ += self.t3cof * t3 + t4 * (self.t4cof + t * self.t5cof)

        nm = self.no
        em = tle.ecco
        if tempa <= 0.0:
            raise PropagationError("satellite has decayed (a < 0)")
        am = ((XKE / nm) ** X2O3) * tempa * tempa
        nm = XKE / (am ** 1.5)
        em -= tempe
        if em >= 1.0 or em < -0.001:
            raise PropagationError("eccentricity out of range")
        em = max(em, 1.0e-6)

        mm += self.no * templ
        xlm = mm + argpm + nodem
        nodem %= TAU
        argpm %= TAU
        xlm %= TAU
        mm = (xlm - argpm - nodem) % TAU

        # long-period periodics
        axnl = em * math.cos(argpm)
        temp = 1.0 / (am * (1.0 - em * em))
        aynl = em * math.sin(argpm) + temp * self.aycof
        xl = mm + argpm + nodem + temp * self.xlcof * axnl

        # Kepler's equation
        u = (xl - nodem) % TAU
        eo1 = u
        tem5 = 9999.9
        for _ in range(10):
            if abs(tem5) < 1.0e-12:
                break
            sineo1 = math.sin(eo1)
            coseo1 = math.cos(eo1)
            tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl
            tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5
            if abs(tem5) >= 0.95:
                tem5 = 0.95 if tem5 > 0.0 else -0.95
            eo1 += tem5
        sineo1 = math.sin(eo1)
        coseo1 = math.cos(eo1)

        ecose = axnl * coseo1 + aynl * sineo1
        esine = axnl * sineo1 - aynl * coseo1
        el2 = axnl * axnl + aynl * aynl
        pl = am * (1.0 - el2)
        if pl < 0.0:
            raise PropagationError("semi-latus rectum < 0")

        rl = am * (1.0 - ecose)
        rdotl = math.sqrt(am) * esine / rl
        rvdotl = math.sqrt(pl) / rl
        betal = math.sqrt(1.0 - el2)
        temp = esine / (1.0 + betal)
        sinu = am / rl * (sineo1 - aynl - axnl * temp)
        cosu = am / rl * (coseo1 - axnl + aynl * temp)
        su = math.atan2(sinu, cosu)
        sin2u = (cosu + cosu) * sinu
        cos2u = 1.0 - 2.0 * sinu * sinu

        temp = 1.0 / pl
        temp1 = 0.5 * J2 * temp
        temp2 = temp1 * temp

        # short-period periodics
        mrt = (rl * (1.0 - 1.5 * temp2 * betal * self.x3thm1)
               + 0.5 * temp1 * self.x1mth2 * cos2u)
        su -= 0.25 * temp2 * self.x7thm1 * sin2u
        xnode = nodem + 1.5 * temp2 * self.cosio * sin2u
        xinc = self.inclo + 1.5 * temp2 * self.cosio * self.sinio * cos2u
        mvt = rdotl - nm * temp1 * self.x1mth2 * sin2u / XKE
        rvdot = rvdotl + nm * temp1 * (self.x1mth2 * cos2u
                                       + 1.5 * self.x3thm1) / XKE

        if mrt < 1.0:
            raise PropagationError("satellite has decayed (below surface)")

        # orientation vectors
        sinsu, cossu = math.sin(su), math.cos(su)
        snod, cnod = math.sin(xnode), math.cos(xnode)
        sini, cosi = math.sin(xinc), math.cos(xinc)
        xmx = -snod * cosi
        xmy = cnod * cosi
        ux = xmx * sinsu + cnod * cossu
        uy = xmy * sinsu + snod * cossu
        uz = sini * sinsu
        vx = xmx * cossu - cnod * sinsu
        vy = xmy * cossu - snod * sinsu
        vz = sini * cossu

        r = (mrt * ux * RE, mrt * uy * RE, mrt * uz * RE)
        v = ((mvt * ux + rvdot * vx) * VKMPERSEC,
             (mvt * uy + rvdot * vy) * VKMPERSEC,
             (mvt * uz + rvdot * vz) * VKMPERSEC)
        return r, v

    def at(self, when: dt.datetime):
        tsince = (when - self.tle.epoch).total_seconds() / 60.0
        return self.propagate(tsince)


# ===========================================================================
# time + reference frames
# ===========================================================================

def julian_day(d: dt.datetime) -> float:
    y, m = d.year, d.month
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    jd = (math.floor(365.25 * (y + 4716)) + math.floor(30.6001 * (m + 1))
          + d.day + b - 1524.5)
    frac = (d.hour + (d.minute + (d.second + d.microsecond / 1e6) / 60.0) / 60.0) / 24.0
    return jd + frac


def gmst(when: dt.datetime) -> float:
    """Greenwich mean sidereal time, radians."""
    tut1 = (julian_day(when) - 2451545.0) / 36525.0
    sec = (-6.2e-6 * tut1 ** 3 + 0.093104 * tut1 ** 2
           + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841)
    return (sec * RAD / 240.0) % TAU


def teme_to_ecef(r, when: dt.datetime):
    """Rotate TEME -> ECEF by GMST.  Polar motion is ignored (metres)."""
    g = gmst(when)
    c, s = math.cos(g), math.sin(g)
    return (c * r[0] + s * r[1], -s * r[0] + c * r[1], r[2])


def ecef_to_geodetic(r):
    """ECEF km -> (lat deg, lon deg, alt km) on WGS84."""
    x, y, z = r
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1.0 - WGS84_E2))
    for _ in range(8):
        sl = math.sin(lat)
        n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sl * sl)
        newlat = math.atan2(z + WGS84_E2 * n * sl, p)
        if abs(newlat - lat) < 1e-13:
            lat = newlat
            break
        lat = newlat
    sl = math.sin(lat)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sl * sl)
    alt = p / math.cos(lat) - n if abs(math.cos(lat)) > 1e-9 else abs(z) - n * (1 - WGS84_E2)
    lon = (lon + math.pi) % TAU - math.pi
    return lat * DEG, lon * DEG, alt


def geodetic_to_ecef(lat_deg, lon_deg, alt_km):
    lat, lon = lat_deg * RAD, lon_deg * RAD
    sl, cl = math.sin(lat), math.cos(lat)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * sl * sl)
    return ((n + alt_km) * cl * math.cos(lon),
            (n + alt_km) * cl * math.sin(lon),
            (n * (1.0 - WGS84_E2) + alt_km) * sl)


class Observer:
    def __init__(self, lat, lon, alt_km=0.0, name="here"):
        self.lat, self.lon, self.alt = lat, lon, alt_km
        self.name = name
        self.ecef = geodetic_to_ecef(lat, lon, alt_km)

    def look(self, sat_ecef):
        """-> (az deg, el deg, range km) in the topocentric SEZ frame."""
        rx = sat_ecef[0] - self.ecef[0]
        ry = sat_ecef[1] - self.ecef[1]
        rz = sat_ecef[2] - self.ecef[2]
        lat, lon = self.lat * RAD, self.lon * RAD
        sl, cl = math.sin(lat), math.cos(lat)
        so, co = math.sin(lon), math.cos(lon)
        south = sl * co * rx + sl * so * ry - cl * rz
        east = -so * rx + co * ry
        zen = cl * co * rx + cl * so * ry + sl * rz
        rng = math.sqrt(rx * rx + ry * ry + rz * rz)
        el = math.asin(max(-1.0, min(1.0, zen / rng))) * DEG
        az = (math.atan2(east, -south) * DEG) % 360.0
        return az, el, rng


# ===========================================================================
# the observation a single satellite produces at a single instant
# ===========================================================================

class Fix:
    __slots__ = ("tle", "az", "el", "rng", "rate", "lat", "lon", "alt",
                 "speed", "error")

    def __init__(self, tle, az=0.0, el=-90.0, rng=0.0, rate=0.0,
                 lat=0.0, lon=0.0, alt=0.0, speed=0.0, error=None):
        self.tle, self.error = tle, error
        self.az, self.el, self.rng, self.rate = az, el, rng, rate
        self.lat, self.lon, self.alt, self.speed = lat, lon, alt, speed

    @property
    def visible(self):
        return self.error is None and self.el > 0.0


def observe(prop: SGP4, obs: Observer, when: dt.datetime) -> Fix:
    try:
        r, v = prop.at(when)
    except PropagationError as exc:
        return Fix(prop.tle, error=str(exc))
    except (ValueError, OverflowError, ZeroDivisionError) as exc:
        return Fix(prop.tle, error=type(exc).__name__)

    ecef = teme_to_ecef(r, when)
    lat, lon, alt = ecef_to_geodetic(ecef)
    az, el, rng = obs.look(ecef)

    # range rate by one-second finite difference, good to ~mm/s here
    try:
        r2, _ = prop.at(when + dt.timedelta(seconds=1))
        _, _, rng2 = obs.look(teme_to_ecef(r2, when + dt.timedelta(seconds=1)))
        rate = rng2 - rng
    except PropagationError:
        rate = 0.0

    speed = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return Fix(prop.tle, az, el, rng, rate, lat, lon, alt, speed)


# ===========================================================================
# pass prediction
# ===========================================================================

def _el_at(prop, obs, when):
    try:
        r, _ = prop.at(when)
    except PropagationError:
        return -90.0
    return obs.look(teme_to_ecef(r, when))[1]


def _refine(prop, obs, lo, hi, rising):
    """Bisect the horizon crossing between lo and hi down to ~1 s."""
    for _ in range(24):
        mid = lo + (hi - lo) / 2
        up = _el_at(prop, obs, mid) > 0.0
        if up == rising:
            hi = mid
        else:
            lo = mid
        if (hi - lo).total_seconds() < 1.0:
            break
    return lo + (hi - lo) / 2


def find_passes(prop: SGP4, obs: Observer, start: dt.datetime,
                hours=48, step=30, min_el=0.0, limit=10):
    """Coarse-step the horizon, then bisect each crossing."""
    passes = []
    t = start
    end = start + dt.timedelta(hours=hours)
    delta = dt.timedelta(seconds=step)
    prev_el = _el_at(prop, obs, t)
    aos = None
    peak_el, peak_t = -90.0, t

    while t < end and len(passes) < limit:
        t += delta
        el = _el_at(prop, obs, t)
        if aos is None and prev_el <= 0.0 < el:
            aos = _refine(prop, obs, t - delta, t, True)
            peak_el, peak_t = el, t
        elif aos is not None:
            if el > peak_el:
                peak_el, peak_t = el, t
            if el <= 0.0 < prev_el:
                los = _refine(prop, obs, t - delta, t, False)
                # sharpen the peak with a short fine sweep
                fine = peak_t - delta
                while fine <= peak_t + delta:
                    e = _el_at(prop, obs, fine)
                    if e > peak_el:
                        peak_el, peak_t = e, fine
                    fine += dt.timedelta(seconds=5)
                if peak_el >= min_el:
                    passes.append({
                        "aos": aos, "los": los, "peak": peak_t,
                        "max_el": peak_el,
                        "aos_az": obs.look(teme_to_ecef(prop.at(aos)[0], aos))[0],
                        "los_az": obs.look(teme_to_ecef(prop.at(los)[0], los))[0],
                        "duration": (los - aos).total_seconds(),
                    })
                aos, peak_el = None, -90.0
        prev_el = el
    return passes


# ===========================================================================
# TLE fetching + cache
# ===========================================================================

def fetch_tles(group="stations", force=False, quiet=False):
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, f"{group}.tle")
    fresh = (os.path.exists(path)
             and time.time() - os.path.getmtime(path) < CACHE_TTL)

    if fresh and not force:
        with open(path, encoding="utf-8") as fh:
            return parse_tle_text(fh.read()), f"cache ({group})"

    url = GP_URL.format(group)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "i_dont_wanna_die/1.0"})
        with urllib.request.urlopen(req, timeout=20) as resp:
            text = resp.read().decode("utf-8", "replace")
        if "1 " not in text:
            raise ValueError("celestrak returned no elements")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        return parse_tle_text(text), f"celestrak ({group})"
    except (urllib.error.URLError, OSError, ValueError) as exc:
        if os.path.exists(path):
            if not quiet:
                print(f"warn: fetch failed ({exc}); using stale cache",
                      file=sys.stderr)
            with open(path, encoding="utf-8") as fh:
                return parse_tle_text(fh.read()), f"stale cache ({group})"
        raise SystemExit(f"could not fetch TLEs and no cache at {path}: {exc}")


def load_config():
    try:
        with open(CONFIG, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def save_config(cfg):
    try:
        os.makedirs(os.path.dirname(CONFIG), exist_ok=True)
        with open(CONFIG, "w", encoding="utf-8") as fh:
            json.dump(cfg, fh, indent=2)
    except OSError:
        pass


# ===========================================================================
# formatting helpers
# ===========================================================================

COMPASS = ("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
           "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")


def compass(az):
    return COMPASS[int((az % 360.0) / 22.5 + 0.5) % 16]


def hhmmss(td_seconds):
    s = int(abs(td_seconds))
    return f"{s // 3600:02d}:{(s % 3600) // 60:02d}:{s % 60:02d}"


def key_for(i):
    return "123456789ABCDEFGHJKLMNPQRSTUVWXYZ"[i] if i < 33 else "*"


# ===========================================================================
# curses UI
# ===========================================================================

def run_tui(stdscr, props, obs, source, group):
    import curses

    curses.curs_set(0)
    stdscr.nodelay(True)
    use_color = curses.has_colors()
    if use_color:
        curses.start_color()
        curses.use_default_colors()
        for n, col in enumerate((curses.COLOR_GREEN, curses.COLOR_YELLOW,
                                 curses.COLOR_BLUE, curses.COLOR_CYAN,
                                 curses.COLOR_RED, curses.COLOR_MAGENTA), 1):
            curses.init_pair(n, col, -1)

    def C(n):
        return curses.color_pair(n) if use_color else 0

    def put(y, x, text, attr=0, width=None, pad=False):
        """Clip to the window and write.  Truncates by character, never by
        byte, so the degree signs cannot smear the layout."""
        h, w = stdscr.getmaxyx()
        if not 0 <= y < h or not 0 <= x < w:
            return
        avail = (w - x) if width is None else min(width, w - x)
        if y == h - 1:
            avail = min(avail, w - x - 1)   # last cell would scroll
        if avail <= 0:
            return
        s = text[:avail]
        if pad:
            s = s.ljust(avail)
        try:
            stdscr.addstr(y, x, s, attr)
        except curses.error:
            pass

    # column budget: head 21 | el/az 19 | range 10 | rate 7 | geo 25
    W_HEAD, W_ELAZ, W_RNG, W_RATE, W_GEO = 21, 19, 10, 7, 25

    def header_text(tw):
        s = f"{'':3}{'SATELLITE':<18}"
        if tw >= W_HEAD + W_ELAZ:
            s += f"{'ELEV':>7}{'AZIM':>8}{'':4}"
        if tw >= W_HEAD + W_ELAZ + W_RNG:
            s += f"{'RANGE':>10}"
        if tw >= W_HEAD + W_ELAZ + W_RNG + W_RATE:
            s += f"{'RATE':>7}"
        if tw >= W_HEAD + W_ELAZ + W_RNG + W_RATE + W_GEO:
            s += f"{'LAT':>8}{'LON':>9}{'ALT':>8}"
        return s

    def row_text(i, f, tw, selected):
        mark = ">" if selected else " "
        s = f"{mark}{key_for(i)} {f.tle.name[:18]:<18}"
        if f.error:
            return s + " " + f.error
        if tw >= W_HEAD + W_ELAZ:
            s += f"{f.el:>6.1f}°{f.az:>7.1f}° {compass(f.az):<3}"
        if tw >= W_HEAD + W_ELAZ + W_RNG:
            s += f"{f.rng:>8.0f}km"
        if tw >= W_HEAD + W_ELAZ + W_RNG + W_RATE:
            s += f"{f.rate:>+7.2f}"
        if tw >= W_HEAD + W_ELAZ + W_RNG + W_RATE + W_GEO:
            s += f"{f.lat:>8.2f}{f.lon:>9.2f}{f.alt:>6.0f}km"
        return s

    def draw_radar(x0, rw, top_y, bot_y, fixes, sel):
        cy = (top_y + bot_y) // 2
        ry = min((bot_y - top_y) // 2, 9)
        rx = min(rw // 2 - 3, ry * 2)   # keep E off the very last cell
        if rx < 5 or ry < 4:
            return
        ccx = x0 + rw // 2
        for yy in range(top_y, bot_y + 1):
            put(yy, x0 - 1, "│", C(4) | curses.A_DIM)
        # one row at a time, so the flat top and bottom of each ellipse
        # stay two clean dots instead of a smear
        for frac in (1.0, 2.0 / 3.0, 1.0 / 3.0):
            fry, frx = frac * ry, frac * rx
            for py in range(int(round(cy - fry)), int(round(cy + fry)) + 1):
                if not top_y <= py <= bot_y:
                    continue
                dy = (cy - py) / fry
                if abs(dy) > 1.0:
                    continue
                dx = frx * math.sqrt(max(0.0, 1.0 - dy * dy))
                put(py, int(round(ccx - dx)), ".", C(4) | curses.A_DIM)
                put(py, int(round(ccx + dx)), ".", C(4) | curses.A_DIM)
        put(cy, ccx, "+", C(4) | curses.A_DIM)
        for lbl, dx, dy in (("N", 0, -1), ("S", 0, 1), ("E", 1, 0), ("W", -1, 0)):
            put(cy + dy * (ry + 1), ccx + dx * (rx + 1), lbl, C(4))
        for i, f in enumerate(fixes):
            if not f.visible:
                continue
            rr = (90.0 - f.el) / 90.0
            a = f.az * RAD
            px = int(round(ccx + rr * rx * math.sin(a)))
            py = int(round(cy - rr * ry * math.cos(a)))
            if top_y <= py <= bot_y:
                if i == sel:
                    attr = C(6) | curses.A_REVERSE | curses.A_BOLD
                elif f.el > 30.0:
                    attr = C(1) | curses.A_BOLD
                else:
                    attr = C(2)
                put(py, px, key_for(i), attr)

    sel_norad = None
    sel = 0
    paused = False
    frozen = None
    radar_on = True
    pass_cache = {}
    HELP = "j/k move  p passes  m radar  r refetch  space pause  q quit"
    status = HELP
    last = 0.0
    fixes = []

    while True:
        now = frozen or dt.datetime.now(dt.timezone.utc)
        if time.time() - last > 0.9 or not fixes:
            pairs = [(observe(p, obs, now), p) for p in props]
            pairs.sort(key=lambda pr: (pr[0].error is not None, -pr[0].el))
            fixes = [p[0] for p in pairs]
            props = [p[1] for p in pairs]
            last = time.time()

        # keep the cursor on the satellite, not on the row number
        if sel_norad is not None:
            for i, f in enumerate(fixes):
                if f.tle.norad == sel_norad:
                    sel = i
                    break
        sel = max(0, min(sel, len(fixes) - 1)) if fixes else 0
        if fixes:
            sel_norad = fixes[sel].tle.norad

        h, w = stdscr.getmaxyx()
        stdscr.erase()

        nvis = sum(1 for f in fixes if f.visible)
        head = (f" i_dont_wanna_die  {now:%Y-%m-%d %H:%M:%S}Z  "
                f"{obs.lat:+.3f},{obs.lon:+.3f} @{obs.alt*1000:.0f}m  "
                f"{nvis}/{len(fixes)} up  [{source}]"
                + ("  PAUSED" if paused else ""))
        put(0, 0, head, C(4) | curses.A_BOLD, pad=True)

        radar_w = 0
        if radar_on and h >= 17 and w - 58 >= 26:
            radar_w = min(34, w - 58)
        table_w = w - radar_w - (1 if radar_w else 0)

        put(2, 0, header_text(table_w), curses.A_UNDERLINE,
            width=table_w, pad=True)

        rows = max(1, h - 6)
        top = max(0, min(sel - rows // 2, len(fixes) - rows))
        for n in range(rows):
            i = top + n
            if i >= len(fixes):
                break
            f = fixes[i]
            if f.error:
                attr = C(5)
            elif f.el > 30.0:
                attr = C(1) | curses.A_BOLD
            elif f.el > 0.0:
                attr = C(2)
            else:
                attr = C(3) | curses.A_DIM
            if i == sel:
                attr |= curses.A_REVERSE
            put(3 + n, 0, row_text(i, f, table_w, i == sel), attr,
                width=table_w, pad=True)

        if radar_w:
            draw_radar(w - radar_w, radar_w, 3, h - 3, fixes, sel)

        if fixes:
            f = fixes[sel]
            t = f.tle
            if f.error:
                info = f"{t.name} [{t.norad}] -- {f.error}"
            else:
                info = (f"{t.name} [{t.norad}]  period {t.period:.1f}min  "
                        f"incl {t.inclo*DEG:.1f}°  "
                        f"TLE age {t.age_days(now):.1f}d  "
                        f"v {f.speed:.2f}km/s")
                nxt = pass_cache.get(t.norad)
                if nxt == "none":
                    info += "  |  no pass in 48h"
                elif nxt:
                    secs = (nxt["aos"] - now).total_seconds()
                    if secs > 0:
                        info += (f"  |  AOS in {hhmmss(secs)} "
                                 f"{compass(nxt['aos_az'])}>"
                                 f"{compass(nxt['los_az'])} "
                                 f"max {nxt['max_el']:.0f}° "
                                 f"for {hhmmss(nxt['duration'])}")
                    else:
                        info += ("  |  OVERHEAD, LOS in "
                                 f"{hhmmss((nxt['los'] - now).total_seconds())}")
            put(h - 2, 0, info, C(6), pad=True)
        put(h - 1, 0, status, C(4) | curses.A_DIM, pad=True)
        stdscr.refresh()

        try:
            ch = stdscr.getch()
        except curses.error:
            ch = -1
        if ch == -1:
            time.sleep(0.1)
            continue

        def move(delta):
            nonlocal sel, sel_norad
            sel = max(0, min(sel + delta, len(fixes) - 1)) if fixes else 0
            sel_norad = fixes[sel].tle.norad if fixes else None

        if ch in (ord("q"), ord("Q"), 27):
            break
        elif ch in (curses.KEY_DOWN, ord("j")):
            move(1)
        elif ch in (curses.KEY_UP, ord("k")):
            move(-1)
        elif ch == curses.KEY_NPAGE:
            move(10)
        elif ch == curses.KEY_PPAGE:
            move(-10)
        elif ch in (curses.KEY_HOME, ord("g")):
            move(-len(fixes))
        elif ch in (curses.KEY_END, ord("G")):
            move(len(fixes))
        elif ch == ord("m"):
            radar_on = not radar_on
        elif ch == ord(" "):
            paused = not paused
            frozen = now if paused else None
        elif ch in (ord("p"), ord("P")) and fixes and not fixes[sel].error:
            t = fixes[sel].tle
            put(h - 1, 0, f"computing passes for {t.name} ...", C(2), pad=True)
            stdscr.refresh()
            ps = find_passes(SGP4(t), obs, now, hours=48, limit=1, min_el=1.0)
            pass_cache[t.norad] = ps[0] if ps else "none"
            status = HELP
        elif ch in (ord("r"), ord("R")):
            put(h - 1, 0, "refetching TLEs ...", C(2), pad=True)
            stdscr.refresh()
            try:
                tles, source = fetch_tles(group, force=True, quiet=True)
                props = [SGP4(t) for t in tles]
                pass_cache.clear()
                last = 0.0
                status = HELP
            except SystemExit as exc:
                status = f"refetch failed: {exc}"


# ===========================================================================
# plain-text output
# ===========================================================================

def print_once(props, obs, now, show_all):
    fixes = sorted((observe(p, obs, now) for p in props),
                   key=lambda f: (f.error is not None, -f.el))
    print(f"{now:%Y-%m-%d %H:%M:%S} UTC   observer {obs.lat:+.4f},{obs.lon:+.4f} "
          f"@ {obs.alt*1000:.0f} m")
    print(f"{'SATELLITE':<24}{'ELEV':>7}{'AZIM':>8}{'':5}{'RANGE':>10}"
          f"{'RATE':>9}{'LAT':>9}{'LON':>10}{'ALT':>8}")
    print("-" * 80)
    shown = 0
    for f in fixes:
        if f.error:
            if show_all:
                print(f"{f.tle.name[:24]:<24}  -- {f.error}")
            continue
        if not show_all and f.el <= 0.0:
            continue
        shown += 1
        print(f"{f.tle.name[:24]:<24}{f.el:>6.1f}°{f.az:>7.1f}° "
              f"{compass(f.az):<4}{f.rng:>9.0f}km{f.rate:>+8.2f}"
              f"{f.lat:>9.2f}{f.lon:>10.2f}{f.alt:>7.0f}km")
    if not shown and not show_all:
        print("(nothing above the horizon -- use --all to see everything)")


def print_passes(props, obs, now, query, hours, min_el):
    q = query.lower()
    hits = [p for p in props if q in p.tle.name.lower() or q == p.tle.norad]
    if not hits:
        print(f"no satellite matching {query!r}")
        return
    for p in hits[:5]:
        print(f"\n{p.tle.name}  [{p.tle.norad}]  "
              f"TLE age {p.tle.age_days(now):.1f} d")
        if p.deep_space:
            print("  deep-space orbit, not supported by this propagator")
            continue
        ps = find_passes(p, obs, now, hours=hours, min_el=min_el, limit=20)
        if not ps:
            print(f"  no pass above {min_el:.0f}° in the next {hours} h")
            continue
        print(f"  {'AOS (UTC)':<21}{'DUR':>8}{'MAX EL':>8}{'FROM':>6}{'TO':>5}")
        for x in ps:
            print(f"  {x['aos']:%Y-%m-%d %H:%M:%S}  {hhmmss(x['duration']):>8}"
                  f"{x['max_el']:>7.1f}°{compass(x['aos_az']):>6}"
                  f"{compass(x['los_az']):>5}")


# ===========================================================================

def main(argv=None):
    cfg = load_config()
    ap = argparse.ArgumentParser(
        prog="i_dont_wanna_die",
        description="zero-dependency satellite tracker (SGP4 from scratch)")
    ap.add_argument("--lat", type=float, default=cfg.get("lat"),
                    help="observer latitude, degrees north")
    ap.add_argument("--lon", type=float, default=cfg.get("lon"),
                    help="observer longitude, degrees east")
    ap.add_argument("--alt", type=float, default=cfg.get("alt", 0.0),
                    help="observer altitude, metres")
    ap.add_argument("--group", default=cfg.get("group", "stations"),
                    help="celestrak group: stations, visual, active, "
                         "starlink, noaa, gps-ops, ...")
    ap.add_argument("--tle", help="read a local TLE file instead of fetching")
    ap.add_argument("--refresh", action="store_true", help="force a TLE refetch")
    ap.add_argument("--once", action="store_true", help="print one table and exit")
    ap.add_argument("--all", action="store_true",
                    help="with --once, include satellites below the horizon")
    ap.add_argument("--passes", metavar="NAME",
                    help="print upcoming passes for a satellite and exit")
    ap.add_argument("--hours", type=float, default=48.0,
                    help="pass-search window (default 48)")
    ap.add_argument("--min-el", type=float, default=10.0,
                    help="minimum peak elevation for --passes (default 10)")
    args = ap.parse_args(argv)

    if args.lat is None or args.lon is None:
        ap.error("need --lat and --lon the first time (they get remembered)")

    obs = Observer(args.lat, args.lon, args.alt / 1000.0)

    if args.tle:
        with open(args.tle, encoding="utf-8") as fh:
            tles = parse_tle_text(fh.read())
        source = os.path.basename(args.tle)
    else:
        tles, source = fetch_tles(args.group, force=args.refresh,
                                  quiet=not (args.once or args.passes))
    if not tles:
        raise SystemExit("no TLEs parsed")

    # only now that the group/file actually worked: a typo'd --group must
    # not get baked into the config and break every later run
    save_config({"lat": args.lat, "lon": args.lon,
                 "alt": args.alt, "group": args.group})

    props = [SGP4(t) for t in tles]
    now = dt.datetime.now(dt.timezone.utc)

    if args.passes:
        print_passes(props, obs, now, args.passes, args.hours, args.min_el)
        return 0
    if args.once:
        print_once(props, obs, now, args.all)
        return 0

    try:
        locale.setlocale(locale.LC_ALL, "")
    except locale.Error:
        pass
    import curses
    try:
        curses.wrapper(run_tui, props, obs, source, args.group)
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
