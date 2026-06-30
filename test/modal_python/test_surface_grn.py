#!/Users/yusong/miniforge3/bin/python

import os
import subprocess
import tempfile
import unittest

import numpy as np
from obspy import read

import pygrt


REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TEST_DIR = os.path.join(REPO_DIR, "test")
MILROW = os.path.join(TEST_DIR, "milrow")
GRT_BIN = os.path.join(REPO_DIR, "pygrt", "C_extension", "bin", "grt")


class TestPythonSurfaceGreenFunctions(unittest.TestCase):
    def setUp(self):
        self.modarr = np.loadtxt(MILROW)

    def test_surface_grn_matches_cli_modsum_sac(self):
        model = pygrt.PyModel1D(self.modarr, 2.0, 1.0)
        disp = model.compute_dispersion(
            freq_range=(0.0, 0.2, 0.05),
            wave_type="Rayleigh",
            nmode=0,
        )

        streams = model.compute_surface_grn(
            disp,
            100.0,
            modes=[0],
            upsampling_n=2,
            calc_upar=True,
        )

        self.assertEqual(len(streams), 1)
        stream = streams[0]
        names = {tr.stats.sac.kcmpnm for tr in stream}
        self.assertIn("EXZ", names)
        self.assertIn("zEXZ", names)
        self.assertIn("rEXZ", names)

        env = os.environ.copy()
        env["PATH"] = os.pathsep.join([os.path.dirname(GRT_BIN), env.get("PATH", "")])
        env["DYLD_LIBRARY_PATH"] = os.pathsep.join(
            ["/Users/yusong/miniforge3/lib", env.get("DYLD_LIBRARY_PATH", "")]
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            phase = os.path.join(tmpdir, "phase_R.nc")
            subprocess.run(
                [
                    GRT_BIN,
                    "eigenv",
                    f"-M{MILROW}",
                    "-F0/0.2/0.05",
                    "-SR",
                    "-N0",
                    f"-C{phase}",
                    "-s",
                ],
                check=True,
                cwd=tmpdir,
                env=env,
            )
            subprocess.run(
                [
                    GRT_BIN,
                    "modsum",
                    f"-C{phase}",
                    "-D2/1",
                    "-R100",
                    "-N0",
                    "-OGRN_NM_0",
                    "-W2",
                    "-e",
                ],
                check=True,
                cwd=tmpdir,
                env=env,
            )
            outdir = os.path.join(tmpdir, "GRN_NM_0", "milrow_2_1_100")

            for kcmpnm in ["EXZ", "EXR", "HFZ", "zEXZ", "rEXZ"]:
                cli_trace = read(os.path.join(outdir, f"{kcmpnm}.sac"))[0]
                py_trace = stream.select(channel=kcmpnm)[0]
                self.assertEqual(py_trace.stats.npts, cli_trace.stats.npts)
                self.assertAlmostEqual(py_trace.stats.delta, cli_trace.stats.delta)
                self.assertEqual(py_trace.stats.sac.kcmpnm, cli_trace.stats.sac.kcmpnm)
                np.testing.assert_allclose(
                    py_trace.data,
                    cli_trace.data,
                    rtol=1e-5,
                    atol=1e-8,
                )

    def test_surface_grn_requires_uniform_nonzero_dispersion(self):
        model = pygrt.PyModel1D(self.modarr, 2.0, 1.0)
        disp = pygrt.PyDispersion(
            "Rayleigh",
            np.array([0.05, 0.12]),
            np.array([[3.0], [3.0]]),
            np.array([[0], [0]]),
            np.array([1, 1]),
        )

        with self.assertRaises(ValueError):
            model.compute_surface_grn(disp, 100.0)


if __name__ == "__main__":
    unittest.main()
