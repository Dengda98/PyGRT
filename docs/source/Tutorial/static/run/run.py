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
# 设置 recv_points 来传入任意点坐标文件
pymod.static_syn(
    recv_points="rcv_pts.txt",
    scale=1e24, strike=33, dip=90, rake=0, output_path="stsyn_points.nc", zne=True,
)
# END SYN POINTS
# ---------------------------------------------------------------------------------

# # ---------------------------------------------------------------------------------
# # 只是备份，可供参考
# import matplotlib.pyplot as plt
# from typing import Union

# def plot_static(static_syn:dict, out:Union[str,None]=None):
#     vars_ = static_syn["variables"]
#     north = vars_["north"]["data"]
#     east = vars_["east"]["data"]
#     Z = vars_["Z"]["data"]
#     N = vars_["N"]["data"]
#     E = vars_["E"]["data"]

#     fig, ax = plt.subplots(1, 1, figsize=(10, 8))
#     # 设计对称色标
#     m = np.max(np.abs(Z)) * 1.2
#     pcm = ax.pcolormesh(east, north, Z, cmap='bwr', vmin=-m, vmax=m)
#     ax.quiver(east, north, E, N, angles='uv', pivot='mid')
#     ax.set_ylim([north[0], north[-1]])
#     ax.set_xlim([east[0], east[-1]])
#     ax.set_aspect('equal')
#     cbar = fig.colorbar(pcm, ax=ax, label='Z(cm)')
#     cbar.formatter.set_powerlimits((0, 0))

#     if out is not None:
#         fig.savefig(out, bbox_inches='tight')
# # ---------------------------------------------------------------------------------


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
