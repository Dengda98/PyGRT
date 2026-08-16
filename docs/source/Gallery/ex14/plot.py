import pygrt
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from obspy import read

# 定义半无限空间模型 
Vp = 8     # km/s
Vs = 4.62  # km/s
Rho = 3.3  # g/cm^3

# 模型数组，半无限空间
modarr = np.array([
    [0.,  Vp, Vs, Rho, 9e10, 9e10],
])

# 剪切模量
mu = Vs**2 * Rho

depsrc = 0.0  # 震源深度 km
deprcv = 0.0  # 台站深度 km

rs = np.array([10]) # 震中距数组，km

nt = 1000     # 总点数，不要求2的幂次
dt = 0.005    # 采样时间间隔(s)

idx = 0      # 震中距索引
r = rs[idx]

t = np.arange(0, nt)*dt * Vs/r 

modfile = "_halfspace_mod"
np.savetxt(modfile, modarr)
pymod = pygrt.PyModel1D(grn="GRN", modelpath=modfile)
# 计算格林函数（仅一个震中距，可用通配符读回）
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, dists=rs, nt=nt, dt=dt,
    converg_method='none',
)
st_none = read("GRN/*/*.sac")
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, dists=rs, nt=nt, dt=dt,
    converg_method='DCM',
)
st_dcm = read("GRN/*/*.sac")
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, dists=rs, nt=nt, dt=dt,
    converg_method='PTAM',
)
st_ptam = read("GRN/*/*.sac")

# 卷积阶跃函数
pygrt.utils.stream_integral(st_none)
pygrt.utils.stream_integral(st_dcm)
pygrt.utils.stream_integral(st_ptam)

coef = np.pi * np.pi * mu * r 

fig, axs = plt.subplots(5, 3, figsize=(12, 10), sharex=True, gridspec_kw=dict(wspace=0.2))
for i, ch in enumerate(['VFR', 'VFZ', 'HFR', 'HFZ', 'HFT']):
    ax = axs[i, 0]
    data = st_none.select(channel=ch)[0].data * coef
    ax.plot(t, data)
    ax.set_ylim(-2, 2)
    ax.set_xlim(0, 2)
    ax.grid()
    ax.text(0.05, 0.92, ch, transform=ax.transAxes, ha='left', va='top', fontsize=12, bbox=dict(facecolor='w'))

    ax = axs[i, 1]
    data = st_dcm.select(channel=ch)[0].data * coef
    ax.plot(t, data)
    ax.set_ylim(-2, 2)
    ax.set_xlim(0, 2)
    ax.grid()
    ax.text(0.05, 0.92, ch, transform=ax.transAxes, ha='left', va='top', fontsize=12, bbox=dict(facecolor='w'))

    ax = axs[i, 2]
    data = st_ptam.select(channel=ch)[0].data * coef
    ax.plot(t, data)
    ax.set_ylim(-2, 2)
    ax.set_xlim(0, 2)
    ax.grid()
    ax.text(0.05, 0.92, ch, transform=ax.transAxes, ha='left', va='top', fontsize=12, bbox=dict(facecolor='w'))

axs[0,0].set_title("None")
axs[0,1].set_title("Apply DCM")
axs[0,2].set_title("Apply PTAM")

fig.savefig("lamb1_compare.svg", bbox_inches='tight')

# 删除中间计算结果，仅保留成图
import shutil
for name in ["GRN", "_halfspace_mod"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
