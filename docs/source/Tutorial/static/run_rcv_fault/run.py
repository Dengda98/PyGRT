import numpy as np
import pygrt

# ------------------------------------------------------------------
# BEGIN GRN
# 指定模型路径以及静态格林函数库的保存路径
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

# 设置多个震源深度
pymod.compute_static_grn(
    depsrc=5.0,
    deprcv=np.arange(0, 10+1e-8, 0.5),
    dists=np.arange(0, 20+1e-8, 0.5),
    calc_upar=True
)
# END GRN
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# BEGIN SYN
# 设置 rcv_fault 来传入 Coulomb 格式的有限断层文件来设置接收点
pymod.compute_static_syn(
    rcv_fault="rcv_fault.inp",
    scale=1e24, strike=33, dip=90, rake=0, output_path="stsyn_rf.nc",
    zne=True, calc_upar=True
)
# END SYN
# ------------------------------------------------------------------
