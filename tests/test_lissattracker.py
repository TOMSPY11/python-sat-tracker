"""Tests for lissattracker.

The SGP4 case is the official verification data from Vallado et al.,
"Revisiting Spacetrack Report #3" (AIAA 2006-6753), satellite 00005.
"""

import datetime as dt
import math
import unittest

import lissattracker as lt


# official verification vectors: tsince (min), position (km), velocity (km/s)
SAT_00005 = (
    "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753",
    "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667",
)
SAT_00005_EXPECTED = [
    (0.0, 7022.46529266, -1400.08296755, 0.03995155,
     1.893841015, 6.405893759, 4.534807250),
    (360.0, -7154.03120202, -3783.17682504, -3536.19412294,
     4.741887409, -4.151817765, -2.093935425),
    (720.0, -7134.59340119, 6531.68641334, 3260.27186483,
     -4.113793027, -2.911922039, -2.557327851),
    (1080.0, 5568.53901181, 4492.06992591, 3863.87641983,
     -4.209106476, 5.159719888, 2.744852980),
    (1440.0, -938.55923943, -6268.18748831, -4294.02924751,
     7.536105209, -0.427127707, 0.989878080),
]


class TestSGP4(unittest.TestCase):
    """The propagator must reproduce the reference vectors to the precision
    they are published at (8 decimal places, i.e. ~10 nm)."""

    def setUp(self):
        self.prop = lt.SGP4(lt.TLE("SAT 00005", *SAT_00005))

    def test_matches_reference_vectors(self):
        for tsince, rx, ry, rz, vx, vy, vz in SAT_00005_EXPECTED:
            with self.subTest(tsince=tsince):
                r, v = self.prop.propagate(tsince)
                self.assertLess(math.dist(r, (rx, ry, rz)), 1e-6,
                                "position drifted from the reference")
                self.assertLess(math.dist(v, (vx, vy, vz)), 1e-7,
                                "velocity drifted from the reference")

    def test_propagation_is_deterministic(self):
        self.assertEqual(self.prop.propagate(123.0), self.prop.propagate(123.0))

    def test_specific_angular_momentum_is_conserved(self):
        """No drag is modelled over a single orbit at this altitude, so
        |r x v| should stay put to a few parts in a thousand."""
        mags = []
        for tsince in range(0, 133, 11):
            r, v = self.prop.propagate(float(tsince))
            h = (r[1] * v[2] - r[2] * v[1],
                 r[2] * v[0] - r[0] * v[2],
                 r[0] * v[1] - r[1] * v[0])
            mags.append(math.dist(h, (0, 0, 0)))
        self.assertLess((max(mags) - min(mags)) / min(mags), 5e-3)


class TestTLEParsing(unittest.TestCase):
    def test_fields(self):
        tle = lt.TLE("SAT 00005", *SAT_00005)
        self.assertEqual(tle.norad, "00005")
        self.assertEqual(tle.name, "SAT 00005")
        self.assertAlmostEqual(tle.inclo * lt.DEG, 34.2682, places=4)
        self.assertAlmostEqual(tle.nodeo * lt.DEG, 348.7242, places=4)
        self.assertAlmostEqual(tle.ecco, 0.1859667, places=7)
        self.assertAlmostEqual(tle.argpo * lt.DEG, 331.7664, places=4)
        self.assertAlmostEqual(tle.mo * lt.DEG, 19.3264, places=4)
        self.assertAlmostEqual(tle.bstar, 0.28098e-4, places=12)

    def test_epoch_is_utc(self):
        tle = lt.TLE("SAT 00005", *SAT_00005)
        self.assertEqual(tle.epoch.tzinfo, dt.timezone.utc)
        self.assertEqual(tle.epoch.year, 2000)

    def test_two_digit_year_pivot(self):
        """57-99 is the 1900s, 00-56 is the 2000s."""
        old = SAT_00005[0][:18] + "57" + SAT_00005[0][20:]
        self.assertEqual(lt.TLE("x", old, SAT_00005[1]).epoch.year, 1957)

    def test_exponential_field(self):
        self.assertAlmostEqual(lt._expfield(" 28098-4"), 0.28098e-4, places=12)
        self.assertAlmostEqual(lt._expfield("-11606-4"), -0.11606e-4, places=12)
        self.assertAlmostEqual(lt._expfield(" 00000-0"), 0.0)
        self.assertAlmostEqual(lt._expfield("        "), 0.0)

    def test_parses_three_line_and_bare_two_line(self):
        three = "ISS\n" + "\n".join(SAT_00005)
        self.assertEqual(lt.parse_tle_text(three)[0].name, "ISS")
        bare = "\n".join(SAT_00005)
        self.assertEqual(lt.parse_tle_text(bare)[0].name, "00005")

    def test_garbage_yields_nothing(self):
        self.assertEqual(lt.parse_tle_text("hello world\nnot a tle"), [])

    def test_deep_space_is_flagged(self):
        """Period >= 225 min is out of scope and must be refused, not guessed."""
        geo = ("1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753",
               "2 00005  34.2682 348.7242 0000667 331.7664  19.3264  1.00270000413667")
        tle = lt.TLE("GEO", *geo)
        self.assertTrue(tle.deep_space)
        with self.assertRaises(lt.PropagationError):
            lt.SGP4(tle).propagate(0.0)


class TestFrames(unittest.TestCase):
    def test_geodetic_round_trip(self):
        for lat, lon, alt in [(0.0, 0.0, 0.0), (51.4779, -0.0015, 0.045),
                              (-33.86, 151.21, 0.058), (89.0, 179.0, 400.0),
                              (-89.0, -179.0, 35786.0)]:
            with self.subTest(lat=lat, lon=lon):
                back = lt.ecef_to_geodetic(lt.geodetic_to_ecef(lat, lon, alt))
                self.assertAlmostEqual(back[0], lat, places=7)
                self.assertAlmostEqual(back[1], lon, places=7)
                self.assertAlmostEqual(back[2], alt, places=6)

    def test_gmst_advances_one_sidereal_day(self):
        """A solar day advances GMST by 360.9856 deg, not 360 -- that extra
        ~0.9856 deg is what makes the stars rise 4 minutes earlier each day.
        Taken mod a full turn, the excess is what should show up."""
        t0 = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)
        excess = math.degrees((lt.gmst(t0 + dt.timedelta(days=1))
                               - lt.gmst(t0)) % lt.TAU)
        self.assertAlmostEqual(excess, 0.9856, places=3)

    def test_gmst_known_epoch(self):
        """J2000.0: GMST is 280.46 deg."""
        j2000 = dt.datetime(2000, 1, 1, 12, 0, 0, tzinfo=dt.timezone.utc)
        self.assertAlmostEqual(math.degrees(lt.gmst(j2000)), 280.46, places=2)

    def test_teme_to_ecef_preserves_radius(self):
        r = (4000.0, -5000.0, 3000.0)
        t = dt.datetime(2026, 6, 1, 7, 30, tzinfo=dt.timezone.utc)
        self.assertAlmostEqual(math.dist(r, (0, 0, 0)),
                               math.dist(lt.teme_to_ecef(r, t), (0, 0, 0)),
                               places=9)


class TestObserver(unittest.TestCase):
    def test_overhead_satellite_is_at_the_zenith(self):
        obs = lt.Observer(13.75, 100.5, 0.0)
        above = lt.geodetic_to_ecef(13.75, 100.5, 500.0)
        az, el, rng = obs.look(above)
        self.assertAlmostEqual(el, 90.0, places=4)
        self.assertAlmostEqual(rng, 500.0, places=3)

    def test_azimuth_points_the_right_way(self):
        obs = lt.Observer(0.0, 0.0, 0.0)
        for d_lat, d_lon, expect in [(10, 0, 0.0), (0, 10, 90.0),
                                     (-10, 0, 180.0), (0, -10, 270.0)]:
            with self.subTest(expect=expect):
                target = lt.geodetic_to_ecef(d_lat, d_lon, 500.0)
                az = obs.look(target)[0]
                self.assertAlmostEqual(az, expect, places=3)

    def test_antipodal_satellite_is_far_below(self):
        obs = lt.Observer(0.0, 0.0, 0.0)
        az, el, rng = obs.look(lt.geodetic_to_ecef(0.0, 180.0, 400.0))
        self.assertLess(el, -80.0)
        self.assertGreater(rng, 12000.0)


class TestHelpers(unittest.TestCase):
    def test_compass(self):
        for az, pt in [(0, "N"), (45, "NE"), (90, "E"), (180, "S"),
                       (270, "W"), (359, "N"), (360, "N"), (22.5, "NNE")]:
            self.assertEqual(lt.compass(az), pt, f"az={az}")

    def test_hhmmss(self):
        self.assertEqual(lt.hhmmss(0), "00:00:00")
        self.assertEqual(lt.hhmmss(61), "00:01:01")
        self.assertEqual(lt.hhmmss(3661), "01:01:01")
        self.assertEqual(lt.hhmmss(-61), "00:01:01")

    def test_keys_are_unique(self):
        keys = [lt.key_for(i) for i in range(33)]
        self.assertEqual(len(set(keys)), 33)


class TestPasses(unittest.TestCase):
    def test_finds_a_pass_and_the_geometry_is_consistent(self):
        """Any pass found must rise, peak and set in that order, with the
        peak elevation at least as high as the endpoints."""
        prop = lt.SGP4(lt.TLE("SAT 00005", *SAT_00005))
        start = prop.tle.epoch
        obs = lt.Observer(34.0, 0.0, 0.0)   # under the inclination
        passes = lt.find_passes(prop, obs, start, hours=24, min_el=0.0, limit=3)
        self.assertGreater(len(passes), 0, "expected at least one pass")
        for p in passes:
            self.assertLess(p["aos"], p["peak"])
            self.assertLess(p["peak"], p["los"])
            self.assertGreater(p["max_el"], 0.0)
            self.assertAlmostEqual(p["duration"],
                                   (p["los"] - p["aos"]).total_seconds(),
                                   places=6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
