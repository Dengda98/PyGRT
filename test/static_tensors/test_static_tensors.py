import numpy as np
import pygrt

depsrc=2
deprcv=0
norths = np.arange(-3., 3.1, 0.2)
easts = np.arange(-2., 2.1, 0.2)

modname="../milrow"

modarr = np.loadtxt(modname)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)
stgrn = pymod.compute_static_grn(norths, easts, calc_upar=True)

stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22, calc_upar=True)

strain = pygrt.utils.compute_strain(stsyn)
stress = pygrt.utils.compute_stress(stsyn)
rotation = pygrt.utils.compute_rotation(stsyn)

stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22, ZNE=True, calc_upar=True)
strain = pygrt.utils.compute_strain(stsyn)
stress = pygrt.utils.compute_stress(stsyn)
rotation = pygrt.utils.compute_rotation(stsyn)