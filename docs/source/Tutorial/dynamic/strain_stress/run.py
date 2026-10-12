# --------------------------------------------------------------------------------------
# BEGIN GRN
import numpy as np
import pygrt
from obspy import read

pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")

# 传入 calc_upar=True 计算空间导数
pymod.greenfn(depsrc=2.0, deprcv=0.0, dists=[30], nt=1500, dt=0.02, calc_upar=True)
stgrn = read("GRN/*/*.sac")
print(stgrn.__str__(extended=True))
# 45 Trace(s) in Stream:
# .SYN..EXZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:29.980000Z | 50.0 Hz, 1500 samples
# ...
# .SYN..zEXZ | ...
# .SYN..rEXZ | ...
# END GRN
# --------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------
# BEGIN SYN DC
# 传入 calc_upar=True 计算空间导数，integrate_order=1 得到阶跃源响应
# time_function 卷积上升 0.5 s、下降 1.5 s 的非对称余弦波
# ?.sac 为位移，[zrt]?.sac 为空间导数
pymod.syn(dist=30.0, azimuth=30.0, scale=1e24, output_path="syn_dc",
          strike=33, dip=50, rake=120, calc_upar=True, integrate_order=1, time_function="c/0.5/1.5")
stsyn = read("syn_dc/?.sac") + read("syn_dc/[zrt]?.sac")
print(stsyn)
# 12 Trace(s) in Stream:
# .SYN..Z  | ...
# .SYN..zZ | ...
# .SYN..rZ | ...
# .SYN..tZ | ...
# END SYN DC
# --------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------
# BEGIN ZNE
# 传入 zne=True 可返回 ZNE 分量，integrate_order=1 得到阶跃源响应
# ?.sac 为位移，[zne]?.sac 为空间导数
pymod.syn(dist=30.0, azimuth=30.0, scale=1e24, output_path="syn_dc_zne",
          strike=33, dip=50, rake=120, calc_upar=True, zne=True, integrate_order=1, time_function="c/0.5/1.5")
stsyn = read("syn_dc_zne/?.sac") + read("syn_dc_zne/[zne]?.sac")
print(stsyn)
# 12 Trace(s) in Stream:
# .SYN..Z  | ...
# .SYN..zZ | ...
# .SYN..nZ | ...
# .SYN..eZ | ...
# END ZNE
# --------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
# BEGIN plot func
from obspy import *
import matplotlib.pyplot as plt

def plot6(st6:Stream, title:str, out:str|None=None, *, zero_channels:tuple[str, ...]=()):
    nt = st6[0].stats.npts
    dt = st6[0].stats.delta
    t = np.arange(nt)*dt

    MAX, MIN = 0, 9e30
    for tr in st6:
        d = np.abs(tr.data)
        if MAX < np.max(d):
            MAX = np.max(d)
        if MIN > np.min(d):
            MIN = np.min(d)

    travtP = st6[0].stats.sac['t0']
    travtS = st6[0].stats.sac['t1']

    fig, axs = plt.subplots(len(st6), 1, figsize=(10, 1.2*len(st6)), gridspec_kw=dict(hspace=0.0), sharex=True)
    for i in range(len(st6)):
        tr = st6[i]
        ax = axs[i]

        ax.plot(t, tr.data, c='k', lw=0.5, label=tr.stats.channel)
        ax.legend(loc='upper left')

        # 理论为零或幅值极小的分量采用整体量级设置纵轴，避免放大数值残差
        m = np.max(np.abs(tr.data))
        ylims = ax.get_ylim()
        if tr.stats.channel in zero_channels or m / MAX < 1e-5:
            ylims = np.array([-1, 1]) * MAX
            ax.set_ylim(ylims)

        # 绘制到时
        ax.vlines(travtP, *ylims, colors='b')
        ax.text(travtP, ylims[1], "P", ha='left', va='top', color='b')
        ax.vlines(travtS, *ylims, colors='r')
        ax.text(travtS, ylims[1], "S", ha='left', va='top', color='r')

        ax.margins(x=0)
        ax.set_xlim([t[0], t[-1]])
        ax.set_ylim(np.array(ylims)*1.2)

    fig.suptitle(title)

    if out is not None:
        fig.savefig(out, bbox_inches='tight')
    return fig, axs
# END plot func
# --------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------
# BEGIN STRAIN
pygrt.utils.strain("syn_dc_zne")
st_strain = read("syn_dc_zne/strain_*.sac")
print(st_strain)
# 6 Trace(s) in Stream:
# .SYN..ZZ | ...
# .SYN..ZN | ...
# .SYN..ZE | ...
# .SYN..NN | ...
# .SYN..NE | ...
# .SYN..EE | ...
# END STRAIN
plot6(st_strain, "Strain", "strain.svg")
# --------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
# BEGIN ROTATION
pygrt.utils.rotation("syn_dc_zne")
st_rotation = read("syn_dc_zne/rotation_*.sac")
print(st_rotation)
# 3 Trace(s) in Stream:
# .SYN..ZN | ...
# .SYN..ZE | ...
# .SYN..NE | ...
# END ROTATION
plot6(st_rotation, "Rotation", "rotation.svg")
# --------------------------------------------------------------------------------------



# BEGIN STRESS
pygrt.utils.stress("syn_dc_zne")
st_stress = read("syn_dc_zne/stress_*.sac")
print(st_stress)
# 6 Trace(s) in Stream:
# .SYN..ZZ | ...
# .SYN..ZN | ...
# .SYN..ZE | ...
# .SYN..NN | ...
# .SYN..NE | ...
# .SYN..EE | ...
# END STRESS
plot6(st_stress, "Stress", "stress.svg", zero_channels=("ZZ", "ZN", "ZE"))
# --------------------------------------------------------------------------------------

# BEGIN SPROJ
# 接收断层形态与震源机制独立设置
pygrt.utils.sproj("syn_dc_zne", strike=90, dip=60, rake=0)
# END SPROJ

# BEGIN COULOMB
pygrt.utils.coulomb("syn_dc_zne", friction=0.4)
# END COULOMB

# BEGIN COULOMB PLOT
st_projected = Stream()
for name, label in zip(("tau_s", "sigma_n", "coulomb"), (r"$\Delta\tau_s$", r"$\Delta\sigma_n$", r"$\Delta CFS$")):
    tr = read(f"syn_dc_zne/{name}.sac")[0]
    tr.data *= 1e-7
    tr.stats.channel = f"{label} (MPa)"
    st_projected += tr
fig, axes = plot6(st_projected, "Stress projection and Coulomb stress change")
axes[-1].set_xlabel("Time (s)")
fig.savefig("coulomb.svg", bbox_inches="tight")
plt.close(fig)
# END COULOMB PLOT

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in ["GRN", "syn_dc", "syn_dc_zne"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
