import numpy as np
from obspy import *
import matplotlib.pyplot as plt
import pygrt

st1 = read("syn/?.sac")
st2 = read("syn_nosoil/?.sac")
nt = st1[0].stats.npts
dt = st1[0].stats.delta
t = np.arange(0, nt)*dt

norm = np.linalg.norm(st1.select(component='R')[0].data, np.inf)


fig, axs = plt.subplots(3, 2, figsize=(12, 6))

for i, (ch, ylabel) in enumerate(zip(['Z', 'R', 'T'], ['Vertical', 'Radial', 'Transverse'])):
    tr = st1.select(component=ch)[0]
    tr2 = st2.select(component=ch)[0]
    ax = axs[i][0]
    ax2 = axs[i][1]
    ax.plot(t, tr.data/norm, 'k-', lw=1.2)
    ax2.plot(t, tr2.data/norm, 'k-', lw=1.2)

    ax.set_xlim([0, 3])
    ax2.set_xlim([0, 3])

    ax.set_ylabel(ylabel, fontsize=13)

    ax.set_ylim([-1, 1])
    ax2.set_ylim([-1, 1])

    ax.grid()
    ax2.grid()


axs[0,0].set_title("with 5m soft surface layer", fontsize=16)
axs[0,1].set_title("without the soft surface layer", fontsize=16)

axs[-1,0].set_xlabel("Time (s)", fontsize=13)
axs[-1,1].set_xlabel("Time (s)", fontsize=13)

fig.savefig('compare.svg', bbox_inches='tight')