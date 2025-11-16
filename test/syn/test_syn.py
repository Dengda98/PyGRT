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

stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22)
stsyn = pygrt.utils.gen_syn_from_gf_SF(stgrn, 1e16, fN=-1, fE=2, fZ=-4, az=22)
stsyn = pygrt.utils.gen_syn_from_gf_DC(stgrn, 1e20, strike=33, dip=44, rake=55, az=22)
stsyn = pygrt.utils.gen_syn_from_gf_MT(stgrn, 1e20, [1, -2, -5, 0.5, 3, 1.2], az=22)


stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22)

sigs = pygrt.sigs.gen_parabola_wave(0.6, dt)
stsyn0 = stsyn.copy()
pygrt.utils.stream_convolve(stsyn0, sigs)
sigs = pygrt.sigs.gen_trap_wave(0.2, 0.4, 0.7, dt)
stsyn0 = stsyn.copy()
pygrt.utils.stream_convolve(stsyn0, sigs)
sigs = pygrt.sigs.gen_triangle_wave(0.8, dt)
stsyn0 = stsyn.copy()
pygrt.utils.stream_convolve(stsyn0, sigs)
sigs = pygrt.sigs.gen_ricker_wave(1.2, dt)
stsyn0 = stsyn.copy()
pygrt.utils.stream_convolve(stsyn0, sigs)
sigs = np.array([0.0, 0.1, 0.2, 0.4, 0.4, 0.4, 0.2, 0.1, 0.0])
stsyn0 = stsyn.copy()
pygrt.utils.stream_convolve(stsyn0, sigs)

stsyn0 = stsyn.copy()
pygrt.utils.stream_integral(stsyn0)
stsyn0 = stsyn.copy()
pygrt.utils.stream_diff(stsyn0)

stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22, calc_upar=True)
stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22, ZNE=True)
stsyn = pygrt.utils.gen_syn_from_gf_EX(stgrn, 1e20, 22, ZNE=True, calc_upar=True)