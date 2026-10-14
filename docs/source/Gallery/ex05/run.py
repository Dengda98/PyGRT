import numpy as np
import matplotlib.pyplot as plt
from obspy import Stream, read
import shutil
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
        fig.tight_layout()
        fig.savefig(out, bbox_inches='tight')
    plt.close(fig)

# 固定同一剪切源，仅改变面积归一化的矩率时间函数，持续时间均为 0.6 s
# integrate_order=1 对应 CLI -I1，得到阶跃型断层错动的位移
moment = 1e24
azimuth = 39.2
time_functions = [
    ("trig",   "t/0.3/0/0.3"),
    ("trap",   "t/0.1/0.2/0.3"),
    ("para",   "p/0.6"),
    ("cosine", "c/0.1/0.5"),
]
for name, time_function in time_functions:
    output_path = f"syn_{name}"
    pymod.syn(dist=dist, azimuth=azimuth, scale=moment, output_path=output_path,
              strike=77, dip=88, rake=99, integrate_order=1, time_function=time_function)
    st = read(f"{output_path}/?.sac")
    sigs = read(f"{output_path}/sig.sac")[0].data
    plot_syn(st, f"{name}.svg", sigs)

# 删除中间计算结果，仅保留成图
shutil.rmtree("GRN")
for name, _ in time_functions:
    shutil.rmtree(f"syn_{name}")
