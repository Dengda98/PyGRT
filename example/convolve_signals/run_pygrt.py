import numpy as np
import pygrt 

dist=10
depsrc=2
deprcv=0.5

nt=1024
dt=0.01

modname="milrow"

modarr = np.loadtxt(modname)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)

# compute green functions
st_grn = pymod.compute_grn(dist, nt, dt)[0]

# synthetic
S=1e24
az=39.2
st = pygrt.utils.gen_syn_from_gf_EXP(st_grn, S, az)
sigs = pygrt.sigs.gen_triangle_wave(0.4, dt)
pygrt.utils.stream_convolve(st, sigs)
pygrt.utils.stream_write_sac(st, "syn_exp/pout")

st = pygrt.utils.gen_syn_from_gf_SF(st_grn, S, 2, -1, 4, az)
sigs = pygrt.sigs.gen_trap_wave(0.1, 0.3, 0.6, dt)
pygrt.utils.stream_convolve(st, sigs)
pygrt.utils.stream_write_sac(st, "syn_sf/pout")

st = pygrt.utils.gen_syn_from_gf_DC(st_grn, S, 77, 88, 99, az)
sigs = pygrt.sigs.gen_parabola_wave(0.6, dt)
pygrt.utils.stream_convolve(st, sigs)
pygrt.utils.stream_write_sac(st, "syn_dc/pout")

st = pygrt.utils.gen_syn_from_gf_MT(st_grn, S, [1,-2,-5,0.5,3,1.2], az)
sigs = pygrt.sigs.gen_ricker_wave(3, dt)
pygrt.utils.stream_convolve(st, sigs)
pygrt.utils.stream_write_sac(st, "syn_mt/pout")