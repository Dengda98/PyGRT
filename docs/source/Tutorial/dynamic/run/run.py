# -----------------------------------------------------------------------------------
# START BUILD MODEL
import numpy as np
import pygrt

# 直接使用模型文件路径
pymod = pygrt.PyModel1D(modelpath="milrow")

# END BUILD MODEL
# -----------------------------------------------------------------------------------



# -----------------------------------------------------------------------------------
# BEGIN GRN
from obspy import read

pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")

# 结果写入 GRN/milrow_{depsrc}_{deprcv}_{dist}/
pymod.compute_grn(
    depsrc=2.0,
    deprcv=0.0,
    dists=[5, 8, 10],
    nt=500,
    dt=0.02,
)
# 需要时再显式读回；多震中距时需指定子目录，单震中距可用 GRN/*/*.sac
stgrn = read("GRN/milrow_2_0_5/*.sac")

print(stgrn)
# 15 Trace(s) in Stream:
# .SYN..EXZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..VFZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# ...
# END GRN
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN plot func
from obspy import *
import matplotlib.pyplot as plt
from typing import Union

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
        fig.savefig(out, bbox_inches='tight')

def plot_int_dif(stsyn:Stream, stsyn_int:Stream, stsyn_dif:Stream, chnl:str, out:Union[str,None]=None):
    nt = stsyn[0].stats.npts
    dt = stsyn[0].stats.delta
    t = np.arange(nt)*dt

    travtP = stsyn[0].stats.sac['t0']
    travtS = stsyn[0].stats.sac['t1']

    fig, axs = plt.subplots(3, 1, figsize=(10, 4), gridspec_kw=dict(hspace=0.0), sharex=True)
    for i, (st, suffix) in enumerate(zip([stsyn, stsyn_int, stsyn_dif], ["", "_int", "_dif"])):
        tr = st.select(channel=chnl)[0]

        ax = axs[i]
        ax.plot(t, tr.data, c='k', lw=0.5, label=f"{tr.stats.channel}{suffix}")
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
            fig.savefig(out, bbox_inches='tight')

# END plot func
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN REUSE GRN
# 若仅使用已算好的格林函数做合成，构造时只需指定 grn，无需 modelpath
pymod = pygrt.PyModel1D(grn="GRN")
# END REUSE GRN
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN SYN EX
# 接之前的代码，使用震中距 10 km 的格林函数
# ?.sac 匹配位移三分量文件名（Z/R/T 或 Z/N/E）
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e24,
    output_path="syn_ex",
    source="EX",
)
stsyn = read("syn_ex/?.sac")
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..Z | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..R | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..T | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_ex.svg")
# END SYN EX
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# BEGIN SYN SF
# 接之前的代码
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e16,
    output_path="syn_sf",
    source="SF",
    force=(1, -0.5, 2),
)
stsyn = read("syn_sf/?.sac")
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..Z | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..R | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..T | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_sf.svg")
# END SYN SF
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN SYN DC
# 接之前的代码
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e24,
    output_path="syn_dc",
    source="DC",
    strike=33,
    dip=50,
    rake=120,
)
stsyn = read("syn_dc/?.sac")
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..Z | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..R | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..T | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_dc.svg")
# END SYN DC
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN SYN TS
# 接之前的代码
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e24,
    output_path="syn_ts",
    source="TS",
    strike=33,
    dip=50,
)
stsyn = read("syn_ts/?.sac")
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..Z | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..R | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..T | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_ts.svg")
# END SYN TS
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN SYN MT
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e24,
    output_path="syn_mt",
    source="MT",
    moment_tensor=(0.1, -0.2, 1.0, 0.3, -0.5, -2.0),
)
stsyn = read("syn_mt/?.sac")
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..Z | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..R | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..T | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_mt.svg")
# END SYN MT
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN ZNE
# 接之前的代码
# 设置 zne=True 可返回 ZNE 分量
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e24,
    output_path="syn_dc_zne",
    source="DC",
    strike=33,
    dip=50,
    rake=120,
    zne=True,
)
stsyn = read("syn_dc_zne/?.sac")
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..Z | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..N | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..E | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_dc_zne.svg")
# END ZNE
# -----------------------------------------------------------------------------------



# -----------------------------------------------------------------------------------
# BEGIN TIME FUNC
# time_function 对应 CLI -D；t1=t2 时梯形波退化为三角波
# 卷积用的时间函数会以 sig.sac 保存在输出目录
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e16,
    output_path="syn_sf_trig",
    source="SF",
    force=(1, -0.5, 2),
    time_function="t/0.3/0.3/0.6",
)
stsyn = read("syn_sf_trig/?.sac")
trig = read("syn_sf_trig/sig.sac")[0].data
plot_syn(stsyn, "syn_sf_trig.svg", trig)
# END TIME FUNC
# -----------------------------------------------------------------------------------



# -----------------------------------------------------------------------------------
# BEGIN INT DIF
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e24,
    output_path="syn_mt_intdif",
    source="MT",
    moment_tensor=(0.1, -0.2, 1.0, 0.3, -0.5, -2.0),
)
stsyn = read("syn_mt_intdif/?.sac")

# 使用 inplace=False，防止原地修改
stsyn_int = pygrt.utils.stream_integral(stsyn, inplace=False)
stsyn_dif = pygrt.utils.stream_diff(stsyn, inplace=False)

for ch in ['Z', 'R', 'T']:
    plot_int_dif(stsyn, stsyn_int, stsyn_dif, ch, f"syn_mt_intdif_{ch}.svg")
# END INT DIF
# -----------------------------------------------------------------------------------

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in [
    "GRN",
    "syn_ex", "syn_sf", "syn_dc", "syn_ts", "syn_mt",
    "syn_dc_zne", "syn_sf_trig", "syn_mt_intdif",
]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
