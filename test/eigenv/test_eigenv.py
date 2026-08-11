#!/Users/yusong/miniforge3/bin/python

import numpy as np

import pygrt


def check_dispersion(wave_type):
    modarr = np.loadtxt("../milrow")
    pymod = pygrt.PyModel1D(modarr, depsrc=2.0, deprcv=3.0)

    disp = pymod.compute_dispersion(
        freq_range=(0.0, 0.15, 0.05),
        wave_type=wave_type,
        nmode=None,
    )

    if disp.freqs.size == 0:
        raise ValueError(f"{wave_type}: empty frequency array.")
    if disp.c.ndim != 2 or disp.c.shape[0] != disp.freqs.size:
        raise ValueError(f"{wave_type}: invalid phase-velocity shape.")
    if disp.ciref.shape != disp.c.shape:
        raise ValueError(f"{wave_type}: invalid reference-index shape.")
    if not np.any(np.isfinite(disp.c)):
        raise ValueError(f"{wave_type}: no finite phase velocities.")

    freqs, mode0 = disp.mode(0)
    if freqs.shape != disp.freqs.shape or mode0.shape != disp.freqs.shape:
        raise ValueError(f"{wave_type}: invalid mode(0) output.")
    if not np.any(np.isfinite(mode0)):
        raise ValueError(f"{wave_type}: mode(0) has no finite values.")


for wave_type in ("Rayleigh", "Love"):
    check_dispersion(wave_type)
