#!/Users/yusong/miniforge3/bin/python

import os
import unittest

import numpy as np

import pygrt


TEST_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MILROW = os.path.join(TEST_DIR, "milrow")


class TestPythonDispersion(unittest.TestCase):
    def setUp(self):
        self.modarr = np.loadtxt(MILROW)

    def test_rayleigh_dispersion_matches_cli_reference_values(self):
        model = pygrt.PyModel1D(self.modarr, 2.0, 1.0)

        disp = model.compute_dispersion(
            freq_range=(0.0, 0.1, 0.05),
            wave_type="Rayleigh",
            nmode=2,
        )

        self.assertIsInstance(disp, pygrt.PyDispersion)
        self.assertEqual(disp.wave_type, "Rayleigh")
        np.testing.assert_allclose(disp.freqs, [0.05, 0.10])
        np.testing.assert_array_equal(disp.modes, [0, 1, 2])
        np.testing.assert_array_equal(disp.cnum, [1, 2])
        self.assertEqual(disp.c.shape, (2, 3))
        self.assertEqual(disp.ciref.shape, (2, 3))
        self.assertTrue(np.isnan(disp.c[0, 1]))
        self.assertEqual(disp.ciref[0, 1], -1)

        self.assertAlmostEqual(disp.c[0, 0], 3.67261968, places=7)
        self.assertAlmostEqual(disp.c[1, 0], 3.30697205, places=7)
        self.assertAlmostEqual(disp.c[1, 1], 4.64060324, places=7)

        freqs, mode0 = disp.mode(0)
        np.testing.assert_allclose(freqs, disp.freqs)
        np.testing.assert_allclose(mode0, disp.c[:, 0])

    def test_love_dispersion_returns_expected_shape(self):
        model = pygrt.PyModel1D(self.modarr, 2.0, 1.0)

        disp = model.compute_dispersion(
            freqs=[0.05, 0.10],
            wave_type="Love",
            nmode=1,
        )

        self.assertEqual(disp.wave_type, "Love")
        np.testing.assert_allclose(disp.freqs, [0.05, 0.10])
        np.testing.assert_array_equal(disp.modes, [0, 1])
        self.assertEqual(disp.c.shape, (2, 2))
        self.assertTrue(np.all(disp.cnum >= 1))

    def test_freq_range_includes_endpoint_despite_float_roundoff(self):
        model = pygrt.PyModel1D(self.modarr, 2.0, 1.0)

        disp = model.compute_dispersion(
            freq_range=(0.0, 0.15, 0.05),
            wave_type="Rayleigh",
            nmode=0,
        )

        np.testing.assert_allclose(disp.freqs, [0.05, 0.10, 0.15])

    def test_dispersion_validates_inputs(self):
        model = pygrt.PyModel1D(self.modarr, 2.0, 1.0)

        with self.assertRaises(ValueError):
            model.compute_dispersion(freqs=[0.1], freq_range=(0, 1, 0.1))
        with self.assertRaises(ValueError):
            model.compute_dispersion(freqs=[-0.1], wave_type="Rayleigh")
        with self.assertRaises(ValueError):
            model.compute_dispersion(freqs=[0.1], wave_type="Body")
        with self.assertRaises(ValueError):
            model.compute_dispersion(freq_range=(0.0, 0.1, -0.05))

    def test_dispersion_rejects_unsupported_boundaries(self):
        model = pygrt.PyModel1D(
            self.modarr,
            2.0,
            1.0,
            topbound="rigid",
            botbound="halfspace",
        )

        with self.assertRaises(NotImplementedError):
            model.compute_dispersion(freqs=[0.1], wave_type="Rayleigh")


if __name__ == "__main__":
    unittest.main()
