# --------------------------------------------------------------------------------------
# BEGIN GRN
import numpy as np
import pygrt
from obspy import read

pymod = pygrt.PyModel1D("milrow")
pymod.set_dynamic_grn_path("GRN")

# 传入 calc_upar=True 计算空间导数
pymod.compute_grn(
    depsrc=2.0,
    deprcv=0.0,
    distarr=[10],
    nt=500,
    dt=0.02,
    calc_upar=True,
)
stgrn = read("GRN/*/*.sac")
print(stgrn.__str__(extended=True))
# 45 Trace(s) in Stream:
# .SYN..EXZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# ...
# .SYN..zEXZ | ...
# .SYN..rEXZ | ...
# END GRN
# --------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------
# BEGIN SYN DC
# 传入 calc_upar=True 计算空间导数
# ?.sac 为位移，[zrt]?.sac 为空间导数
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e24,
    output_path="syn_dc",
    source="DC",
    strike=33,
    dip=50,
    rake=120,
    calc_upar=True,
)
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
# 传入 zne=True 可返回 ZNE 分量
# ?.sac 为位移，[zne]?.sac 为空间导数
pymod.compute_syn(
    dist=10.0,
    azimuth=30.0,
    scale=1e24,
    output_path="syn_dc_zne",
    source="DC",
    strike=33,
    dip=50,
    rake=120,
    calc_upar=True,
    zne=True,
)
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

def plot6(st6:Stream, title:str, out:str|None=None):
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

        # 相对全局最大值过小时，固定到同样数量级下的坐标轴便于展示近零分量
        m = np.max(np.abs(tr.data))
        ylims = ax.get_ylim()
        if m / MAX < 1e-5:
            ylims = np.array([-1, 1]) * MAX
            ax.set_ylim(ylims)

        # 绘制到时
        ax.vlines(travtP, *ylims, colors='b')
        ax.text(travtP, ylims[1], "P", ha='left', va='top', color='b')
        ax.vlines(travtS, *ylims, colors='r')
        ax.text(travtS, ylims[1], "S", ha='left', va='top', color='r')

        ax.set_xlim([t[0], t[-1]])
        ax.set_ylim(np.array(ylims)*1.2)

    fig.suptitle(title)

    if out is not None:
        fig.savefig(out, bbox_inches='tight')
# END plot func
# --------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------
# BEGIN STRAIN
# return_result=True 时按文件名前缀读回 strain_*.sac，避免与位移等混在一起
st_strain = pygrt.utils.compute_strain("syn_dc_zne", return_result=True)
print(st_strain)
# 6 Trace(s) in Stream:
# .SYN..ZZ | ...
# .SYN..ZN | ...
# .SYN..ZE | ...
# .SYN..NN | ...
# .SYN..NE | ...
# .SYN..EE | ...
plot6(st_strain, "Strain", "strain.svg")
# END STRAIN
# --------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
# BEGIN ROTATION
st_rotation = pygrt.utils.compute_rotation("syn_dc_zne", return_result=True)
print(st_rotation)
# 3 Trace(s) in Stream:
# .SYN..ZN | ...
# .SYN..ZE | ...
# .SYN..NE | ...
plot6(st_rotation, "Rotation", "rotation.svg")
# END ROTATION
# --------------------------------------------------------------------------------------



# BEGIN STRESS
st_stress = pygrt.utils.compute_stress("syn_dc_zne", return_result=True)
print(st_stress)
# 6 Trace(s) in Stream:
# .SYN..ZZ | ...
# .SYN..ZN | ...
# .SYN..ZE | ...
# .SYN..NN | ...
# .SYN..NE | ...
# .SYN..EE | ...
plot6(st_stress, "Stress", "stress.svg")
# END STRESS
# --------------------------------------------------------------------------------------

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in ["GRN", "syn_dc", "syn_dc_zne"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
