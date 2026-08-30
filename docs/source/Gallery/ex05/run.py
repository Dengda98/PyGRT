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

pymod = pygrt.PyModel1D(grn="GRN", modelpath=modname)

# compute green functions
pymod.greenfn(depsrc=depsrc, deprcv=deprcv, dists=[dist], nt=nt, dt=dt)

def plot_syn(stsyn:Stream, out:Union[str,None]=None, sigs:Union[np.ndarray,None]=None):
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
        ax.plot(t[:len(sigs)], sigs, 'k-', lw=0.5)
        axs = axs[1:]

    travtP = traces[0].stats.sac['t0']
    travtS = traces[0].stats.sac['t1']

    for i, tr in enumerate(traces):
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
# ?.sac 匹配位移三分量文件名（Z/R/T）
# time_function 对应 CLI -D
# 所有时间函数使用面积归一化（除雷克子波使用最大幅值为1）
# 卷积用的时间函数保存在输出目录的 sig.sac
S=1e24
az=39.2
pymod.syn(dist=dist, azimuth=az, scale=S, output_path="syn_ex", time_function="t/0.2/0.2/0.4")
st = read("syn_ex/?.sac")
sigs = read("syn_ex/sig.sac")[0].data
plot_syn(st, "trig.svg", sigs)

pymod.syn(dist=dist, azimuth=az, scale=S, output_path="syn_sf", force=(2, -1, 4), time_function="t/0.1/0.3/0.6")
st = read("syn_sf/?.sac")
sigs = read("syn_sf/sig.sac")[0].data
plot_syn(st, "trap.svg", sigs)

pymod.syn(dist=dist, azimuth=az, scale=S, output_path="syn_dc", strike=77, dip=88, rake=99, time_function="p/0.6")
st = read("syn_dc/?.sac")
sigs = read("syn_dc/sig.sac")[0].data
plot_syn(st, "para.svg", sigs)

pymod.syn(dist=dist, azimuth=az, scale=S, output_path="syn_mt", moment_tensor=(1, -2, -5, 0.5, 3, 1.2), time_function="r/3")
st = read("syn_mt/?.sac")
sigs = read("syn_mt/sig.sac")[0].data
plot_syn(st, "rick.svg", sigs)

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in ["GRN", "syn_ex", "syn_sf", "syn_dc", "syn_mt"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
