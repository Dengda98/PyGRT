#!/Users/yusong/miniforge3/bin/python

import numpy as np

import pygrt


modarr = np.loadtxt("../milrow")
pymod = pygrt.PyModel1D(modarr, depsrc=2.0, deprcv=1.0)

disp = pymod.compute_dispersion(
    freq_range=(0.0, 0.2, 0.05),
    wave_type="Rayleigh",
    nmode=0,
)
streams = pymod.compute_surface_grn(
    disp,
    distarr=100.0,
    modes=[0],
    upsampling_n=2,
    calc_upar=True,
)

if len(streams) != 1:
    raise ValueError(f"Expected one stream, got {len(streams)}.")

stream = streams[0]
channel_names = {tr.stats.sac.kcmpnm for tr in stream}
for channel in ("EXZ", "EXR", "zEXZ", "rEXZ"):
    if channel not in channel_names:
        raise ValueError(f"Missing component {channel}.")

    trace = stream.select(channel=channel)[0]
    if trace.stats.npts <= 0:
        raise ValueError(f"{channel}: invalid npts.")
    if trace.stats.delta <= 0.0:
        raise ValueError(f"{channel}: invalid delta.")
    if not np.all(np.isfinite(trace.data)):
        raise ValueError(f"{channel}: non-finite data.")
