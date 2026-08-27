import numpy as np
import pygrt

# 计算格林函数库
pymod = pygrt.PyModel1D(modelpath="prem.flat.20", stgrn="stgrn.nc")
pymod.static_greenfn(
    depsrc=np.arange(1.0, 36.0+1e-8, 2.0),
    deprcv=0.0,
    dists=np.arange(0.0, 500.0+1e-8, 2.0),
    calc_upar=True
)
