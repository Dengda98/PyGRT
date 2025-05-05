# -----------------------------------------------------------------
# BEGIN GRN
import numpy as np
import matplotlib.pyplot as plt
from typing import Union
import pygrt 

modarr = np.loadtxt("mod1")

pymod = pygrt.PyModel1D(modarr, depsrc=0.01, deprcv=0.0)

_ = pymod.compute_grn(distarr=[1], nt=500, dt=0.02, vmin_ref=0.1, Length=20, statsfile="pygrtstats")

vels = np.arange(0.1, 0.6, 0.001)

kerDct = pygrt.utils.read_kernels_freqs("pygrtstats", vels)


Type = "HF_v"
vels = kerDct['_vels']
freqs = kerDct['_freqs']
data = kerDct[Type].copy()
data[...] = data/np.max(np.abs(data), axis=1)[:,None]
fig, ax = plt.subplots(1, 1)
pcm = ax.pcolormesh(freqs, vels, np.abs(np.imag(data)).T, vmin=0, vmax=1, shading='nearest')
fig.colorbar(pcm)

ax.set_xlim([0.5, 25])
ax.set_ylim([0.1, 0.6])
ax.set_xlabel("Frequency (km)")
ax.set_ylabel("Velocity (m/s)")
ax.set_title(f"Im[{Type}]")

fig.tight_layout()
fig.savefig("imag_G.png", dpi=100)