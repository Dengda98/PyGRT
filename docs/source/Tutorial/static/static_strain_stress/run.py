import numpy as np
import pygrt

# BEGIN
# 重新计算格林函数，因为以下计算需要位移空间导数
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

# 计算格林函数
# 传入 calc_upar=True 可计算位移空间导数
pymod.static_greenfn(depsrc=2.0, deprcv=0.0, dists=np.arange(0.0, 10.0 + 1e-8, 0.15), calc_upar=True)

norths = [-3.0, 3.0, 0.15]
easts = [-2.5, 2.5, 0.15]

# 传入 calc_upar=True 可合成位移空间导数
# 传入 zne=True 返回 ZNE 分量
pymod.static_syn(
    norths=norths, easts=easts,
    scale=1e24, strike=33, dip=50, rake=120, output_path="stsyn_dc_zne.nc",
    zne=True, calc_upar=True
)

# 计算应变 / 旋转 / 应力，结果写回同一 nc 文件
pygrt.utils.static_strain("stsyn_dc_zne.nc")
pygrt.utils.static_rotation("stsyn_dc_zne.nc")
pygrt.utils.static_stress("stsyn_dc_zne.nc")
# END

# BEGIN SPROJ
# 接收断层形态与震源机制独立设置
pygrt.utils.static_sproj("stsyn_dc_zne.nc", strike=90, dip=60, rake=0)
# END SPROJ

# BEGIN COULOMB
pygrt.utils.static_coulomb("stsyn_dc_zne.nc", friction=0.4)
# END COULOMB
