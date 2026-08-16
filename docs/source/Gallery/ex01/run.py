import numpy as np
from obspy import *
import matplotlib.pyplot as plt
from typing import Union
import pygrt

pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")
pymod.compute_grn(depsrc=5.0, deprcv=0.0, dists=[180], nt=1400, dt=0.1)
# ?.sac 匹配位移三分量文件名（Z/R/T）
# integrate_order=1 对应 CLI -I1，得到阶跃型位移
pymod.compute_syn(
    dist=180.0, azimuth=39.2, scale=1e24, output_path="syn_dc",
    source="DC", strike=77, dip=88, rake=99,
    integrate_order=1,
)
st = read("syn_dc/?.sac")

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
        tr = stsyn.select(channel=comp)[0]
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

fig, axs = plot_syn(st)
fig.savefig('cover.svg', bbox_inches='tight')

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in ["GRN", "syn_dc"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
