from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pygrt
from matplotlib.colors import Normalize

from plot_fault import read_faults

# ------------------------------------------------------------------
# BEGIN GRN
# 指定模型路径以及静态格林函数库的保存路径
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

# 设置多个震源深度
pymod.compute_static_grn(
    depsrc=np.arange(0, 8+1e-8, 0.5),
    deprcv=0.0,
    dists=np.arange(0, 60+1e-8, 0.5)
)
# END GRN
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# BEGIN SYN
# 设置 src_fault 来传入 Coulomb 格式的有限断层文件
pymod.compute_static_syn(
    norths=(-20, 20, 1),
    easts=(-20, 20, 1),
    src_fault="fault.inp", output_path="stsyn_ff.nc", zne=True,
)
# END SYN
