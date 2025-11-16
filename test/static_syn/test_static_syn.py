import numpy as np
import pygrt

depsrc=2
deprcv=0
xarr = np.arange(-3., 3.1, 0.2)
yarr = np.arange(-2., 2.1, 0.2)

modname="../milrow"

modarr = np.loadtxt(modname)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)
stgrn = pymod.compute_static_grn(xarr, yarr, calc_upar=True)

stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22)
stsyn = pygrt.utils.gen_syn_from_gf_SF(stgrn, 1e16, fN=-1, fE=2, fZ=-4, az=22)
stsyn = pygrt.utils.gen_syn_from_gf_DC(stgrn, 1e20, strike=33, dip=44, rake=55, az=22)
stsyn = pygrt.utils.gen_syn_from_gf_MT(stgrn, 1e20, [1, -2, -5, 0.5, 3, 1.2], az=22)


