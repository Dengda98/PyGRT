# BEGIN GRN
import numpy as np
import pygrt

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

pymod.static_greenfn(depsrc=2.0, deprcv=0.0, dists=np.arange(0.0, 10.0 + 1e-8, 0.1))
# END GRN
