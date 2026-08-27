import numpy as np
import pygrt
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from obspy import read

modname = sys.argv[1]

deprcv = 0.0
depsrc = 5.0
dists = np.arange(0.1, 50, 1)

pymod = pygrt.PyModel1D(grn="GRN", stgrn="stgrn.nc", modelpath=modname)

# 零频频谱: nt=1 时只算 ω=0
# C 反变换写 SAC 时乘了 df=1/(nt*dt)，故时域首点 * dt 即还原频域幅值
dt = 100.0
pymod.greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists.tolist(), nt=1, dt=dt, zeta=1.0, keepAllFreq=True)
# 多震中距：按子目录名末段解析 dist，再按 dists 顺序读入
dist2st = {
    float(p.name.rsplit("_", 1)[-1]): read(str(p / "*.sac"))
    for p in Path("GRN").iterdir() if p.is_dir()
}
stgrnLst = [dist2st[float(d)] for d in dists]

# 静态解
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists.tolist())
static_grn = pygrt.utils.read_static_nc("stgrn.nc")

# 绘制
coef = 1e-20 * 1e25  # 1e25 为地震矩

# 绘制零频结果
fig, ax = plt.subplots(figsize=(6,6))
dynamic_Z = np.zeros_like(dists)
dynamic_R = np.zeros_like(dists)
for i, st in enumerate(stgrnLst):
    # SAC 中 Z 已取反为向上为正，无需再乘 -1
    dynamic_Z[i] = st.select(channel='EXZ')[0].data[0] * dt * coef
    dynamic_R[i] = st.select(channel='EXR')[0].data[0] * dt * coef

ax.plot(dists, dynamic_Z, 'k', label='Dynamic Z')
ax.plot(dists, dynamic_R, 'k', label='Dynamic R')

ax.plot(dists, static_grn['variables']['EXZ']['data'][0, 0, 0] * coef, 'ro', ms=4, label='Static Z')
ax.plot(dists, static_grn['variables']['EXR']['data'][0, 0, 0] * coef, 'bo', ms=4, label='Static R')

ax.set_xlim(0, 50)
ax.set_xlabel('Distance (km)')
ax.set_ylim(ymin=0)
ax.set_ylabel('Displacement (cm)')

ax.legend(loc='upper right')
ax.grid()

fig.savefig(f'{modname}.svg', bbox_inches='tight')

# 删除中间计算结果，仅保留成图
import shutil
for name in ["GRN", "stgrn.nc"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
