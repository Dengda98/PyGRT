import numpy as np
import pygrt 

modarr = np.loadtxt("milrow")

pymod = pygrt.PyModel1D(modarr, 5, 0)

st_grn = pymod.compute_grn(
    distarr=[2500], 
    nt=2000, 
    dt=1, 
    Length=20, 
    safilonTol=1e-2,  # 自适应采样精度
    filonCut=10,
    delayT0=100,
)[0]
