# BEGIN GRN
import numpy as np
import pygrt

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

pymod.static_greenfn(depsrc=2.0, deprcv=0.0, dists=np.arange(0.0, 10.0 + 1e-8, 0.1))
# END GRN

# BEGIN REUSE STGRN
# 若仅使用已算好的静态格林函数做合成，构造时只需指定 stgrn，无需 modelpath
pymod = pygrt.PyModel1D(stgrn="stgrn.nc")
# END REUSE STGRN

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

# BEGIN SYN POINTS
# 设置 rcv_points 来传入任意点坐标文件
pymod.static_syn(
    rcv_points="rcv_pts.txt",
    scale=1e24, strike=33, dip=90, rake=0, output_path="stsyn_points.nc", zne=True,
)
# END SYN POINTS
