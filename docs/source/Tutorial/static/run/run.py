# ---------------------------------------------------------------------------------
# BEGIN GRN
import numpy as np
import pygrt

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

pymod.static_greenfn(depsrc=2.0, deprcv=0.0, dists=np.arange(0.0, 10.0 + 1e-8, 0.1))
# END GRN
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN REUSE STGRN
# 若仅使用已算好的静态格林函数做合成，构造时只需指定 stgrn，无需 modelpath
pymod = pygrt.PyModel1D(stgrn="stgrn.nc")
# END REUSE STGRN
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN GRID
norths = [-3.0, 3.0, 0.15]
easts = [-2.5, 2.5, 0.15]

# 爆炸源
pymod.static_syn(
    norths=norths, easts=easts,
    scale=1e24, output_path="stsyn_ex.nc", zne=True
)

# 单力源
pymod.static_syn(
    norths=norths, easts=easts,
    scale=1e16, force=(1, -0.5, 2), output_path="stsyn_sf.nc", zne=True
)

# 剪切源
pymod.static_syn(
    norths=norths, easts=easts,
    scale=1e24, strike=33, dip=50, rake=120, output_path="stsyn_dc.nc", zne=True
)

# 剪切源
pymod.static_syn(
    norths=norths, easts=easts,
    scale=1e24, strike=33, dip=90, rake=0, output_path="stsyn_dc2.nc", zne=True
)

# 张裂源
pymod.static_syn(
    norths=norths, easts=easts,
    scale=1e24, strike=33, dip=50, output_path="stsyn_ts.nc", zne=True
)

# 张裂源
pymod.static_syn(
    norths=norths, easts=easts,
    scale=1e24, strike=33, dip=90, output_path="stsyn_ts2.nc", zne=True
)

# 矩张量源
pymod.static_syn(
    norths=norths, easts=easts,
    scale=1e24, moment_tensor=(0.1, -0.2, 1.0, 0.3, -0.5, -2.0), output_path="stsyn_mt.nc", zne=True,
)

# 矩张量源
pymod.static_syn(
    norths=norths, easts=easts,
    scale=1e24, moment_tensor=(0, -0.2, 0, 0, 0, 0), output_path="stsyn_mt2.nc", zne=True,
)
# END SYN GRID
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN POINTS
# 设置 rcv_points 来传入任意点坐标文件
pymod.static_syn(
    rcv_points="rcv_pts.txt",
    scale=1e24, strike=33, dip=90, rake=0, output_path="stsyn_points.nc", zne=True,
)
# END SYN POINTS
# ---------------------------------------------------------------------------------

import matplotlib.pyplot as plt

static_syn = pygrt.utils.read_nc_variables("stsyn_points.nc")
fig, ax = pygrt.utils.plot_static_displacement(static_syn, title="Static displacement for -M33/90/0")
ax.plot(0.0, 0.0, marker="*", markersize=12, lw=0, color="black", label="Epicenter")
limit = 1.15 * max(np.max(np.abs(static_syn["north"])), np.max(np.abs(static_syn["east"])))
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.legend(loc="upper right")
fig.savefig("syn_points.svg", dpi=200, bbox_inches="tight")
plt.close(fig)


# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in [
    "stgrn.nc",
    "stsyn_ex.nc", "stsyn_sf.nc", "stsyn_dc.nc", "stsyn_dc2.nc",
    "stsyn_ts.nc", "stsyn_ts2.nc", "stsyn_mt.nc", "stsyn_mt2.nc",
    "stsyn_points.nc",
]:
    p = Path(name)
    if p.is_file():
        p.unlink(missing_ok=True)
