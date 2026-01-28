import numpy as np
import pygrt 
from obspy import *
from compare_func import compare3, update_dict, static_compare3


dist=10
depsrc=2
deprcv=3.3

nt=1024
dt=0.01

modname="../milrow"

modarr = np.loadtxt(modname)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)

AVGRERR = []

#-------------------------- Dynamic -----------------------------------------
# compute green functions
st_grn = pymod.compute_grn(dist, nt, dt, calc_upar=True)[0]

S=1e24
az=39.2

fn=2
fe=-1
fz=4

stk=77
dip=88
rak=99

M11=1
M12=-2
M13=-5
M22=0.5
M23=3
M33=1.2

for ZNE in [False, True]:
    suffix = "-N" if ZNE else ""
    # synthetic
    st = pygrt.utils.gen_syn_from_gf_EX(st_grn, S, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_triangle_wave(0.4, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, f"syn_ex{suffix}/", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    strotation = pygrt.utils.compute_rotation(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, f"syn_ex{suffix}/strain_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(strotation, f"syn_ex{suffix}/rotation_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, f"syn_ex{suffix}/stress_", ZNE=ZNE, dim2=True))
    

    st = pygrt.utils.gen_syn_from_gf_SF(st_grn, S, fn, fe, fz, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_trap_wave(0.1, 0.3, 0.6, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, f"syn_sf{suffix}/", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    strotation = pygrt.utils.compute_rotation(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, f"syn_sf{suffix}/strain_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(strotation, f"syn_sf{suffix}/rotation_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, f"syn_sf{suffix}/stress_", ZNE=ZNE, dim2=True))

    st = pygrt.utils.gen_syn_from_gf_DC(st_grn, S, stk, dip, rak, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_parabola_wave(0.6, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, f"syn_dc{suffix}/", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    strotation = pygrt.utils.compute_rotation(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, f"syn_dc{suffix}/strain_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(strotation, f"syn_dc{suffix}/rotation_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, f"syn_dc{suffix}/stress_", ZNE=ZNE, dim2=True))

    st = pygrt.utils.gen_syn_from_gf_TS(st_grn, S, stk, dip, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_parabola_wave(0.6, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, f"syn_ts{suffix}/", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    strotation = pygrt.utils.compute_rotation(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, f"syn_ts{suffix}/strain_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(strotation, f"syn_ts{suffix}/rotation_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, f"syn_ts{suffix}/stress_", ZNE=ZNE, dim2=True))

    st = pygrt.utils.gen_syn_from_gf_MT(st_grn, S, [M11,M12,M13,M22,M23,M33], az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_ricker_wave(3, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, f"syn_mt{suffix}/", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    strotation = pygrt.utils.compute_rotation(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, f"syn_mt{suffix}/strain_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(strotation, f"syn_mt{suffix}/rotation_", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, f"syn_mt{suffix}/stress_", ZNE=ZNE, dim2=True))


#-------------------------- Static -----------------------------------------
# 为了方便测试，避免引入其他因素的误差，这里有意避开 0
xarr = np.arange(-3.1, 3.2, 0.6)
yarr = np.arange(-4.1, 4.2, 0.8)

static_grn = pymod.compute_static_grn(xarr, yarr, calc_upar=True)
AVGRERR2 = []
# plot_static(static_grn, "static/stgrn.nc")

for ZNE in [False, True]:
    suffix = "-N" if ZNE else ""
    static_syn = pygrt.utils.gen_syn_from_gf_EX(static_grn, S, ZNE=ZNE, calc_upar=True)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_ex{suffix}.nc"))

    static_syn = pygrt.utils.gen_syn_from_gf_SF(static_grn, S, fn, fe, fz, ZNE=ZNE, calc_upar=True)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_sf{suffix}.nc"))

    static_syn = pygrt.utils.gen_syn_from_gf_DC(static_grn, S, stk, dip, rak, ZNE=ZNE, calc_upar=True)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_dc{suffix}.nc"))

    static_syn = pygrt.utils.gen_syn_from_gf_TS(static_grn, S, stk, dip, ZNE=ZNE, calc_upar=True)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_ts{suffix}.nc"))

    static_syn = pygrt.utils.gen_syn_from_gf_MT(static_grn, S, [M11,M12,M13,M22,M23,M33], ZNE=ZNE, calc_upar=True)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_mt{suffix}.nc"))


print("---------------- dynamic --------------------")
AVGRERR = np.array(AVGRERR)
print(AVGRERR)
print(np.mean(AVGRERR), np.min(AVGRERR), np.max(AVGRERR))
print("---------------- static --------------------")
AVGRERR2 = np.array(AVGRERR2)
print(AVGRERR2)
print(np.mean(AVGRERR2), np.min(AVGRERR2), np.max(AVGRERR2))

if np.mean(AVGRERR) > 0.05 or np.mean(AVGRERR2) > 1e-5:
    raise ValueError

