import numpy as np
from obspy import Stream, read
import matplotlib.pyplot as plt

def plot_syn(stsyn: Stream, out: str | None = None, sigs: np.ndarray | None = None):
    traces = list(stsyn)
    order = {ch: i for i, ch in enumerate("ZRTNE")}
    traces.sort(key=lambda tr: order.get(tr.stats.channel, 99))

    figsize = (10, 4)
    nrow = len(traces)
    if sigs is not None:
        nrow += 1
        figsize = (10, 4.5)

    fig, axs = plt.subplots(nrow, 1, figsize=figsize, gridspec_kw=dict(hspace=0.0), sharex=True)
    if nrow == 1:
        axs = [axs]
    nt = traces[0].stats.npts
    dt = traces[0].stats.delta
    t = np.arange(nt)*dt

    if sigs is not None:
        ax = axs[0]
        ax.plot(t[:len(sigs)], sigs, 'k-', lw=0.8)
        axs = axs[1:]

    travtP = traces[0].stats.sac['t0']
    travtS = traces[0].stats.sac['t1']

    for i, tr in enumerate(traces):
        ax = axs[i]
        ax.plot(t, tr.data, c='k', lw=0.8, label=tr.stats.channel)
        ax.legend(loc='upper left')

        ylims = ax.get_ylim()
        # 绘制到时
        ax.vlines(travtP, *ylims, colors='b', linewidths=0.5)
        ax.text(travtP, ylims[1], "P", ha='left', va='top', color='b')
        ax.vlines(travtS, *ylims, colors='r', linewidths=0.5)
        ax.text(travtS, ylims[1], "S", ha='left', va='top', color='r')

        ax.set_xlim([t[0], t[-1]])
        ax.set_ylim(np.array(ylims)*1.2)

    if out is not None:
        fig.savefig(out, bbox_inches='tight')
    plt.close(fig)


def plot_syn_sources(stsyns, labels, out=None):
    fig, axs = plt.subplots(len(stsyns), 1, figsize=(10, 1.4 * len(stsyns)), gridspec_kw=dict(hspace=0.0), sharex=True)
    if len(stsyns) == 1:
        axs = [axs]

    for ax, stsyn, label in zip(axs, stsyns, labels):
        tr = stsyn.select(channel="Z")[0]
        nt = tr.stats.npts
        dt = tr.stats.delta
        t = np.arange(nt) * dt

        ax.plot(t, tr.data, c='k', lw=0.8)
        ylims = ax.get_ylim()

        travtP = tr.stats.sac['t0']
        travtS = tr.stats.sac['t1']
        ax.vlines(travtP, *ylims, colors='b', linewidths=0.5)
        ax.text(travtP, ylims[1], "P", ha='left', va='top', color='b')
        ax.vlines(travtS, *ylims, colors='r', linewidths=0.5)
        ax.text(travtS, ylims[1], "S", ha='left', va='top', color='r')

        ax.text(0.01, 0.85, label, transform=ax.transAxes, ha='left', va='top')
        ax.set_xlim([t[0], t[-1]])
        ax.set_ylim(np.array(ylims) * 1.2)

    axs[-1].set_xlabel("Time (s)")
    if out is not None:
        fig.savefig(out, bbox_inches='tight')
        plt.close(fig)

source_streams = [
    read("syn_ex/?.sac"), read("syn_sf/?.sac"),
    read("syn_dc/?.sac"), read("syn_dc2/?.sac"),
    read("syn_ts/?.sac"), read("syn_mt/?.sac"),
]
source_labels = [
    "-S1e24 -I1", "-F1/-0.5/2", "-M33/50/120 -I1", "-M33/90/0 -I1",
    "-M33/50 -I1", "-T0.1/-0.2/1.0/0.3/-0.5/-2.0 -I1",
]
plot_syn_sources(source_streams, source_labels, "syn_sources.svg")


stsyn = read("syn_dc_zne/?.sac")
plot_syn(stsyn, "syn_dc_zne.svg")

stsyn = read("syn_dc_trig/?.sac")
trig = read("syn_dc_trig/sig.sac")[0].data
plot_syn(stsyn, "syn_dc_trig.svg", trig)

stsyn = read("syn_dc_cosine/?.sac")
cosine = read("syn_dc_cosine/sig.sac")[0].data
plot_syn(stsyn, "syn_dc_cosine.svg", cosine)
