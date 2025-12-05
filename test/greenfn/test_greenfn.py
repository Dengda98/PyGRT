import numpy as np
import pygrt

dist=10
depsrc=2
deprcv=3

nt=600
dt=0.02

modname="../milrow"

modarr = np.loadtxt(modname)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)


stgrn = pymod.compute_grn(dist, nt, dt)[0]
stgrn = pymod.compute_grn(dist, nt, dt, calc_upar=True)[0]
stgrn = pymod.compute_grn(dist, nt, dt, zeta=0.6, upsampling_n=10)[0]
stgrn = pymod.compute_grn(dist, nt, dt, freqband=[1, 10])[0]
stgrn = pymod.compute_grn(dist, nt, dt, Length=20)[0]

stgrn = pymod.compute_grn(dist, nt, dt, k0=4, ampk=1.2, keps=1e-3, vmin_ref=1.5)[0]

stgrn = pymod.compute_grn(2000, 1400, 1.0, Length=20, safilonTol=1e-3)[0]
stgrn = pymod.compute_grn(2000, 1400, 1.0, Length=20, safilonTol=1e-3, keepAllFreq=True)[0]

stgrn = pymod.compute_grn(dist, nt, dt, Length=20, statsfile="GRN_grtstats")[0]
stgrn = pymod.compute_grn(dist, nt, dt, Length=20, statsfile="GRN_grtstats", statsidxs=[1,10,20])[0]

# multi distances
stgrn = pymod.compute_grn([6, 8, 10], nt, dt)[0]
