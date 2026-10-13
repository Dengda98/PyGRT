import numpy as np
import pygrt

# 限定本例所需的深度和震中距范围，使用自适应 Filon 积分加快计算
pymod = pygrt.PyModel1D(modelpath="usgs.model", grn="GRN")
pymod.greenfn(
    depsrc=np.arange(2.0, 35.0 + 1e-8, 3.0),
    deprcv=0.0,
    dists=np.arange(180.0, 300.0 + 1e-8, 10.0),
    nt=500, dt=1.0, upsampling_n=10, keepAllFreq=True,
    safilonTol=1e-2, gf_source=["DC"], calc_upar=True
)
