import numpy as np
import pygrt

# ------------------------------------------------------------------
# BEGIN GRN
# 指定模型路径以及静态格林函数库的保存路径
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

# 设置多个震源深度
pymod.compute_static_grn(
    depsrc=np.arange(0, 8+1e-8, 0.5),
    deprcv=0.0,
    dists=np.arange(0, 60+1e-8, 0.5),
    calc_upar=True
)
# END GRN
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# BEGIN SYN
# 设置 src_fault 来传入 Coulomb 格式的有限断层文件
pymod.compute_static_syn(
    norths=(-20, 20, 1),
    easts=(-20, 20, 1),
    src_fault="fault.inp", output_path="stsyn_ff.nc",
    zne=True, calc_upar=True
)
# END SYN
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# BEGIN COULOMB
# 计算应力张量
pygrt.utils.compute_stress("stsyn_ff.nc")
# 将应力张量投影到指定形态的断层面上
pygrt.utils.compute_sproj("stsyn_ff.nc", strike=59, dip=90, rake=180)
# 指定等效摩擦系数，计算库伦应力
pygrt.utils.compute_coulomb("stsyn_ff.nc", friction=0.75)
# END COULOMB
# ------------------------------------------------------------------
