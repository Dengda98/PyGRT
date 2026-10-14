import numpy as np
import pygrt

# 使用 USGS 配套 FSP 中的六层模型计算静态格林函数库
pymod = pygrt.PyModel1D(modelpath="usgs.model", stgrn="stgrn.nc")
pymod.static_greenfn(
    depsrc=np.arange(1.0, 35.0 + 1e-8, 2.0),
    deprcv=0.0,
    dists=np.arange(0.0, 500.0 + 1e-8, 2.0),
    calc_upar=True
)
