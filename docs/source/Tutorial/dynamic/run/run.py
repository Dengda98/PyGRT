# -----------------------------------------------------------------------------------
# START BUILD MODEL
import numpy as np
import pygrt 

# option 1:
# modarr = np.loadtxt("milrow")

# option 2
modarr = np.array([
    [0.2, 3.4, 1.7, 2.3, 9e10, 9e10],
    [0.6, 3.7, 1.9, 2.4, 9e10, 9e10],
    [0.5, 4.2, 2.1, 2.4, 9e10, 9e10],
    [0.5, 4.6, 2.3, 2.5, 9e10, 9e10],
    [0.7, 4.9, 2.8, 2.6, 9e10, 9e10],
    [0.5, 5.1, 2.9, 2.7, 9e10, 9e10],
    [6.0, 5.9, 3.3, 2.7, 9e10, 9e10],
    [28., 6.9, 4.0, 2.8, 9e10, 9e10],
    [0.,  8.2, 4.7, 3.2, 9e10, 9e10],
])
# END BUILD MODEL
# -----------------------------------------------------------------------------------



# -----------------------------------------------------------------------------------
# BEGIN GRN
modarr = np.loadtxt("milrow")

pymod = pygrt.PyModel1D(modarr, depsrc=2.0, deprcv=0.0)

# 多个震中距的格林函数以列表形式返回，其中每个元素为 |Stream| 类。
stgrnLst = pymod.compute_grn(
    distarr=[5,8,10],
    nt=500, dt=0.02
)

print(stgrnLst[0])
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
        fig.savefig(out, dpi=100)

def plot_int_dif(stsyn:Stream, stsyn_int:Stream, stsyn_dif:Stream, comp:str, out:Union[str,None]=None):
    nt = stsyn[0].stats.npts
    dt = stsyn[0].stats.delta
    t = np.arange(nt)*dt

    travtP = stsyn[0].stats.sac['t0']
    travtS = stsyn[0].stats.sac['t1']

    fig, axs = plt.subplots(3, 1, figsize=(10, 4), gridspec_kw=dict(hspace=0.0), sharex=True)
    for i, (st, suffix) in enumerate(zip([stsyn, stsyn_int, stsyn_dif], ["", "_int", "_dif"])):
        tr = st.select(component=comp)[0]

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
            fig.savefig(out, dpi=100)

# END plot func
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN SYN EXP
# 接之前的代码
idx = 2
stgrn = stgrnLst[idx]   # 选择格林函数

stsyn = pygrt.utils.gen_syn_from_gf_EXP(stgrn, M0=1e24, az=30)
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..EXZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..EXR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..EXT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_exp.png")
# END SYN EXP
# -----------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------
# BEGIN SYN SF
# 接之前的代码
idx = 2
stgrn = stgrnLst[idx]   # 选择格林函数

stsyn = pygrt.utils.gen_syn_from_gf_SF(stgrn, S=1e16, fN=1, fE=-0.5, fZ=2, az=30)
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..SFZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..SFR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..SFT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_sf.png")
# END SYN SF
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN SYN DC
# 接之前的代码
idx = 2
stgrn = stgrnLst[idx]   # 选择格林函数

stsyn = pygrt.utils.gen_syn_from_gf_DC(stgrn, M0=1e24, strike=33, dip=50, rake=120, az=30)
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..DCZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DCR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DCT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_dc.png")
# END SYN DC
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN SYN MT
idx = 2
stgrn = stgrnLst[idx]   # 选择格林函数

stsyn = pygrt.utils.gen_syn_from_gf_MT(stgrn, M0=1e24, MT=[0.1,-0.2,1.0,0.3,-0.5,-2.0], az=30)
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..MTZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..MTR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..MTT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_mt.png")
# END SYN MT
# -----------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# BEGIN ZNE
# 接之前的代码
idx = 2
stgrn = stgrnLst[idx]   # 选择格林函数

# 设置ZNE=True可返回ZNE分量
stsyn = pygrt.utils.gen_syn_from_gf_DC(stgrn, M0=1e24, strike=33, dip=50, rake=120, az=30, ZNE=True)
print(stsyn)
# 3 Trace(s) in Stream:
# .SYN..DCZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DCN | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DCE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
plot_syn(stsyn, "syn_dc_zne.png")
# END ZNE
# -----------------------------------------------------------------------------------



# -----------------------------------------------------------------------------------
# BEGIN TIME FUNC
idx = 2
stgrn = stgrnLst[idx]   # 选择格林函数

stsyn = pygrt.utils.gen_syn_from_gf_SF(stgrn, S=1e16, fN=1, fE=-0.5, fZ=2, az=30)
# 生成时间函数
trig = pygrt.sigs.gen_triangle_wave(0.6, 0.02)
# 卷积，原地修改
pygrt.utils.stream_convolve(stsyn, trig)
plot_syn(stsyn, "syn_sf_trig.png", trig)
# END TIME FUNC
# -----------------------------------------------------------------------------------



# -----------------------------------------------------------------------------------
# BEGIN INT DIF
idx = 2
stgrn = stgrnLst[idx]   # 选择格林函数

stsyn = pygrt.utils.gen_syn_from_gf_MT(stgrn, M0=1e24, MT=[0.1,-0.2,1.0,0.3,-0.5,-2.0], az=30)

# 使用inplace=False，防止原地修改
stsyn_int = pygrt.utils.stream_integral(stsyn, inplace=False)
stsyn_dif = pygrt.utils.stream_diff(stsyn, inplace=False)

for ch in ['Z', 'R', 'T']:
    plot_int_dif(stsyn, stsyn_int, stsyn_dif, ch, f"syn_mt_intdif_{ch}.png")
# END INT DIF
# -----------------------------------------------------------------------------------
