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
stgrn = pymod.compute_grn(dist, nt, dt, calc_upar=True)[0]

stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22, calc_upar=True)

strain = pygrt.utils.compute_strain(stsyn)
stress = pygrt.utils.compute_stress(stsyn)
rotation = pygrt.utils.compute_rotation(stsyn)

stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22, ZNE=True, calc_upar=True)
strain = pygrt.utils.compute_strain(stsyn)
stress = pygrt.utils.compute_stress(stsyn)
rotation = pygrt.utils.compute_rotation(stsyn)


