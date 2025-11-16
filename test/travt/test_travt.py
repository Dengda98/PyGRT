import numpy as np
import pygrt

depsrc=2
deprcv=0

modname="../milrow"

modarr = np.loadtxt(modname)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)

for dist in [2,3,4,5]:
    tp, ts = pymod.compute_travt1d(dist)
    print(dist, tp, ts)