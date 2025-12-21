import numpy as np
import matplotlib.pyplot as plt
from obspy import *
from typing import Union
import pygrt 

dist=10
depsrc=2
deprcv=0.5

nt=1024
dt=0.01

modname="milrow"

modarr = np.loadtxt(modname)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)

# compute green functions
st_grn = pymod.compute_grn(dist, nt, dt)[0]

def plot_syn(stsyn:Stream, out:Union[str,None]=None, sigs:Union[np.ndarray,None]=None):
    figsize = (10, 4)
    nrow = 3
    if sigs is not None:
        nrow += 1
        figsize = (10, 4.5)

    fig, axs = plt.subplots(nrow, 1, figsize=figsize, gridspec_kw=dict(hspace=0.0), sharex=True)
    nt = stsyn[0].stats.npts
    dt = stsyn[0].stats.delta
    t = np.arange(nt)*dt

    if sigs is not None:
        ax = axs[0]
        ax.plot(t[:len(sigs)], sigs, 'k-', lw=0.5)
        axs = axs[1:]

    travtP = stsyn[0].stats.sac['t0']
    travtS = stsyn[0].stats.sac['t1']

    for i, tr in enumerate(stsyn):
        ax = axs[i]
        ax.plot(t, tr.data, c='k', lw=0.5, label=tr.stats.channel)
        ax.legend(loc='upper left')

        ylims = ax.get_ylim()
        # 绘制到时
        ax.vlines(travtP, *ylims, colors='b')
        ax.text(travtP, ylims[1], "P", ha='left', va='top', color='b')
        ax.vlines(travtS, *ylims, colors='r')
        ax.text(travtS, ylims[1], "S", ha='left', va='top', color='r')

        ax.set_xlim([t[0], t[-1]])
        ax.set_ylim(np.array(ylims)*1.2)

    if out is not None:
        fig.tight_layout()
        fig.savefig(out, bbox_inches='tight')

# synthetic
S=1e24
az=39.2
st = pygrt.utils.gen_syn_from_gf_EX(st_grn, S, az)
sigs = pygrt.sigs.gen_triangle_wave(0.4, dt)
pygrt.utils.stream_convolve(st, sigs)
plot_syn(st, "trig.svg", sigs)

st = pygrt.utils.gen_syn_from_gf_SF(st_grn, S, 2, -1, 4, az)
sigs = pygrt.sigs.gen_trap_wave(0.1, 0.3, 0.6, dt)
pygrt.utils.stream_convolve(st, sigs)
plot_syn(st, "trap.svg", sigs)

st = pygrt.utils.gen_syn_from_gf_DC(st_grn, S, 77, 88, 99, az)
sigs = pygrt.sigs.gen_parabola_wave(0.6, dt)
pygrt.utils.stream_convolve(st, sigs)
plot_syn(st, "para.svg", sigs)

st = pygrt.utils.gen_syn_from_gf_MT(st_grn, S, [1,-2,-5,0.5,3,1.2], az)
sigs = pygrt.sigs.gen_ricker_wave(3, dt)
pygrt.utils.stream_convolve(st, sigs)
plot_syn(st, "rick.svg", sigs)
