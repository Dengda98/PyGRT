import numpy as np
from obspy import *
import matplotlib.pyplot as plt
from typing import Union
import pygrt

modarr = np.loadtxt("milrow")
pymod = pygrt.PyModel1D(modarr, 5, 0)
stgrn = pymod.compute_grn(distarr=[180], nt=1400, dt=0.1)[0]
st = pygrt.utils.gen_syn_from_gf_DC(stgrn, M0=1e24, strike=77, dip=88, rake=99, az=39.2)

def plot_syn(stsyn:Stream, sigs:Union[np.ndarray,None]=None):
    figsize = (10, 5.5)
    nrow = 3
    if sigs is not None:
        nrow += 1
        figsize = (10, 4.5)

    fig, axs = plt.subplots(nrow, 1, figsize=figsize, gridspec_kw=dict(hspace=0.0), sharex=True,)
    nt = stsyn[0].stats.npts
    dt = stsyn[0].stats.delta
    t = np.arange(nt)*dt

    if sigs is not None:
        ax = axs[0]
        ax.plot(t[:len(sigs)], sigs, 'k-', lw=0.5)
        axs = axs[1:]

    travtP = stsyn[0].stats.sac['t0']
    travtS = stsyn[0].stats.sac['t1']

    for i, comp in enumerate(['R', 'T', 'Z']):
        ax = axs[i]
        tr = st.select(component=comp)[0]
        ax.plot(t, tr.data, c='k', lw=0.5, label=tr.stats.channel[-1])
        ax.legend(loc='upper left')

        ylims = ax.get_ylim()
        # 绘制到时
        ax.vlines(travtP, *ylims, colors='b')
        ax.text(travtP, ylims[1], "P", ha='left', va='top', color='b')
        ax.vlines(travtS, *ylims, colors='r')
        ax.text(travtS, ylims[1], "S", ha='left', va='top', color='r')

        ax.set_xlim([t[0], t[-1]])
        ax.set_ylim(np.array(ylims)*1.2)
        ax.grid()

    return fig, axs

pygrt.utils.stream_integral(st)
fig, axs = plot_syn(st)
fig.savefig('cover.svg', bbox_inches='tight')