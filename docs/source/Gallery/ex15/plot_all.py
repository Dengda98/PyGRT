import numpy as np
import pygrt
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from obspy import read

modname = sys.argv[1]

deprcv = 0.0
depsrc = 5.0
dists = np.arange(0.1, 50, 0.5)

pymod = pygrt.PyModel1D(grn="GRN", stgrn="stgrn.nc", modelpath=modname)

# 零频频谱: nt=1 时只算 ω=0
# C 反变换写 SAC 时乘了 df=1/(nt*dt)，故时域首点 * dt 即还原频域幅值
dt = 500.0
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, dists=dists.tolist(),
    nt=1, dt=dt, zeta=1.0, keepAllFreq=True,
)
# 多震中距：按子目录名末段解析 dist，再按 dists 顺序读入
dist2st = {
    float(p.name.rsplit("_", 1)[-1]): read(str(p / "*.sac"))
    for p in Path("GRN").iterdir() if p.is_dir()
}
stgrnLst = [dist2st[float(d)] for d in dists]

# 静态解
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, dists=dists.tolist(),
)
static_grn = pygrt.utils.read_static_nc("stgrn.nc")

# 绘制零频结果
fig, axs = plt.subplots(2, 3, figsize=(12, 8), gridspec_kw=dict(hspace=0.3, wspace=0.3)) # 
axs = axs.ravel()

srctypes = ['EX', 'VF', 'HF', 'DD', 'DS', 'SS']
for isrc, (src, src2) in enumerate(zip(srctypes,
                    ['Explosion', 'Vertical Force', 'Horizontal Force', '45°-Dip Slip', '90°-Dip Slip', 'Strike Slip'])):
    ax = axs[isrc]

    if src in ['VF', 'HF']:
        coef = 1e-15 * 1e20  # 1e20 dyne
    else:
        coef = 1e-20 * 1e25  # 1e25 dyne·cm

    dynamic_Z = np.zeros_like(dists)
    dynamic_R = np.zeros_like(dists)
    dynamic_T = np.zeros_like(dists)
    for i, st in enumerate(stgrnLst):
        # SAC 中 Z 已取反为向上为正，无需再乘 -1
        dynamic_Z[i] = st.select(channel=f'{src}Z')[0].data[0] * dt * coef
        dynamic_R[i] = st.select(channel=f'{src}R')[0].data[0] * dt * coef
        if src not in ['EX', 'VF', 'DD']:
            dynamic_T[i] = st.select(channel=f'{src}T')[0].data[0] * dt * coef

    ms = 2
    lw = 0.8
    ax.plot(dists, dynamic_Z, 'k', lw=lw, label='Dynamic Z')
    ax.plot(dists, static_grn['variables'][f'{src}Z']['data'][0, 0, 0] * coef, 'ro', ms=ms, label='Static Z')
    ax.plot(dists, dynamic_R, 'k', lw=lw, label='Dynamic R')
    ax.plot(dists, static_grn['variables'][f'{src}R']['data'][0, 0, 0] * coef, 'bo', ms=ms, label='Static R')
    if src not in ['EX', 'VF', 'DD']:
        ax.plot(dists, dynamic_T, 'k', lw=lw, label='Dynamic T')
        ax.plot(dists, static_grn['variables'][f'{src}T']['data'][0, 0, 0] * coef, 'go', ms=ms, label='Static T')

    ax.set_xlim(0, 50)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Displacement (cm)')
    ax.set_title(src2)
    ax.grid()

ax = axs[2]
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.92), 
          ncols=3, fontsize='large', columnspacing=4.0,
          bbox_transform=fig.transFigure)


fig.savefig(f'{modname}_all.svg', bbox_inches='tight')

# 删除中间计算结果，仅保留成图
import shutil
for name in ["GRN", "stgrn.nc"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
