import numpy as np
import pygrt 
from obspy import *
from compare_func import update_dict, static_compare3

depsrc=2
deprcv=3.3

modname="../milrow"

modarr = np.loadtxt(modname)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)

#-------------------------- Static -----------------------------------------
# 为了方便测试，避免引入其他因素的误差，这里有意避开 0
xarr = np.arange(-3.1, 3.2, 0.6)
yarr = np.arange(-4.1, 4.2, 0.8)
S=1e24

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

static_grn = pymod.compute_static_grn(distarr=np.arange(0, 10+1e-8, 0.1), calc_upar=True)
AVGRERR2 = []
# plot_static(static_grn, "static/stgrn.nc")

for ZNE in [False, True]:
    suffix = "-N" if ZNE else ""
    static_syn = pygrt.utils.gen_syn_from_gf_EX(static_grn, S, ZNE=ZNE, calc_upar=True, xarr=xarr, yarr=yarr)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_ex{suffix}.nc"))

    static_syn = pygrt.utils.gen_syn_from_gf_SF(static_grn, S, fn, fe, fz, ZNE=ZNE, calc_upar=True, xarr=xarr, yarr=yarr)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_sf{suffix}.nc"))

    static_syn = pygrt.utils.gen_syn_from_gf_DC(static_grn, S, stk, dip, rak, ZNE=ZNE, calc_upar=True, xarr=xarr, yarr=yarr)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_dc{suffix}.nc"))

    static_syn = pygrt.utils.gen_syn_from_gf_TS(static_grn, S, stk, dip, ZNE=ZNE, calc_upar=True, xarr=xarr, yarr=yarr)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_ts{suffix}.nc"))

    static_syn = pygrt.utils.gen_syn_from_gf_MT(static_grn, S, [M11,M12,M13,M22,M23,M33], ZNE=ZNE, calc_upar=True, xarr=xarr, yarr=yarr)
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    update_dict(static_syn, ststrain, "strain_")
    update_dict(static_syn, ststress, "stress_")
    update_dict(static_syn, strotation, "rotation_")
    AVGRERR2.append(static_compare3(static_syn, f"static/stsyn_mt{suffix}.nc"))


print("---------------- static --------------------")
AVGRERR2 = np.array(AVGRERR2)
print(AVGRERR2)
print(np.mean(AVGRERR2), np.min(AVGRERR2), np.max(AVGRERR2))

if np.mean(AVGRERR2) > 1e-5:
    raise ValueError
