import numpy as np
import pygrt 
from obspy import *


def compare3(st_py:Stream, c_prefix:str, ZNE:bool=False, dim2:bool=False):
    """return average relative error"""
    
    if ZNE:
        pattern = "[ZNE]"
    else:
        pattern = "[ZRT]"

    if dim2:
        pattern = pattern*2

    st_c = read(f"{c_prefix}{pattern}.sac")

    error = 0.0
    nerr = 0
    for tr_c in st_c:
        tr_py = st_py.select(channel=tr_c.stats.channel)[0]
        rerr = np.mean(np.abs(tr_c.data - tr_py.data) / np.abs(tr_c.data))
        if np.isnan(rerr):
            rerr = 0.0
        error += rerr
        nerr += 1 

    return error/nerr


def static_compare3(resDct:dict, c_prefix:str, ZNE:bool=False, dim2:bool=False):
    """return average relative error"""
    midname = c_prefix.split("_")[-1][:2].upper()

    if ZNE:
        chs2 = ['Z', 'N', 'E']
    else:
        chs2 = ['Z', 'R', 'T']

    chs1 = chs2.copy()

    nx = len(resDct['_xarr'])
    ny = len(resDct['_yarr'])

    c_res = np.loadtxt(c_prefix)
    c_resDct = {}
    if not dim2:
        for i, c in enumerate(chs2):
            c_resDct[f"{midname}{c}"] = c_res[:,i+2].reshape((nx, ny), order='F')
    else:
        chs1 = []
        for i1, c1 in enumerate(chs2):
            for i2, c2 in enumerate(chs2):
                if i2 < i1:
                    continue
                chs1.append(f"{c1}{c2}")

        for i, cc in enumerate(chs1):
            c_resDct[cc] = c_res[:,i+2].reshape((nx, ny), order='F')

    error = 0.0
    nerr = 0

    for k in c_resDct.keys():
        val1 = resDct[k]
        val2 = c_resDct[k]
        rerr = np.mean(np.abs(val1 - val2) / np.abs(val1))
        if np.isnan(rerr):
            rerr = 0.0
        error += rerr
        nerr += 1 

    return error/nerr



dist=10
depsrc=2
deprcv=3.3

nt=1024
dt=0.01

modname="milrow"

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
    # synthetic
    
    st = pygrt.utils.gen_syn_from_gf_EXP(st_grn, S, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_triangle_wave(0.4, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, "syn_exp/cout", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, "syn_exp/cout.strain.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, "syn_exp/cout.stress.", ZNE=ZNE, dim2=True))
    

    st = pygrt.utils.gen_syn_from_gf_SF(st_grn, S, fn, fe, fz, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_trap_wave(0.1, 0.3, 0.6, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, "syn_sf/cout", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, "syn_sf/cout.strain.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, "syn_sf/cout.stress.", ZNE=ZNE, dim2=True))

    st = pygrt.utils.gen_syn_from_gf_DC(st_grn, S, stk, dip, rak, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_parabola_wave(0.6, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, "syn_dc/cout", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, "syn_dc/cout.strain.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, "syn_dc/cout.stress.", ZNE=ZNE, dim2=True))

    st = pygrt.utils.gen_syn_from_gf_MT(st_grn, S, [M11,M12,M13,M22,M23,M33], az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_ricker_wave(3, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, "syn_mt/cout", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, "syn_mt/cout.strain.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, "syn_mt/cout.stress.", ZNE=ZNE, dim2=True))


#-------------------------- Static -----------------------------------------
xarr = np.linspace(-3, 3, 11)
yarr = np.linspace(-3, 3, 11)
static_grn = pymod.compute_static_grn(xarr, yarr, calc_upar=True)

for ZNE in [False, True]:
    suffix = "-N" if ZNE else ""
    static_syn = pygrt.utils.gen_syn_from_gf_EXP(static_grn, S, ZNE=ZNE, calc_upar=True)
    AVGRERR.append(static_compare3(static_syn, f"static/stsyn_exp{suffix}", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    AVGRERR.append(static_compare3(ststrain, f"static/strain_exp{suffix}", ZNE=ZNE, dim2=True))
    AVGRERR.append(static_compare3(ststress, f"static/stress_exp{suffix}", ZNE=ZNE, dim2=True))

    static_syn = pygrt.utils.gen_syn_from_gf_SF(static_grn, S, fn, fe, fz, ZNE=ZNE, calc_upar=True)
    AVGRERR.append(static_compare3(static_syn, f"static/stsyn_sf{suffix}", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    AVGRERR.append(static_compare3(ststrain, f"static/strain_sf{suffix}", ZNE=ZNE, dim2=True))
    AVGRERR.append(static_compare3(ststress, f"static/stress_sf{suffix}", ZNE=ZNE, dim2=True))

    static_syn = pygrt.utils.gen_syn_from_gf_DC(static_grn, S, stk, dip, rak, ZNE=ZNE, calc_upar=True)
    AVGRERR.append(static_compare3(static_syn, f"static/stsyn_dc{suffix}", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    AVGRERR.append(static_compare3(ststrain, f"static/strain_dc{suffix}", ZNE=ZNE, dim2=True))
    AVGRERR.append(static_compare3(ststress, f"static/stress_dc{suffix}", ZNE=ZNE, dim2=True))

    static_syn = pygrt.utils.gen_syn_from_gf_MT(static_grn, S, [M11,M12,M13,M22,M23,M33], ZNE=ZNE, calc_upar=True)
    AVGRERR.append(static_compare3(static_syn, f"static/stsyn_mt{suffix}", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    AVGRERR.append(static_compare3(ststrain, f"static/strain_mt{suffix}", ZNE=ZNE, dim2=True))
    AVGRERR.append(static_compare3(ststress, f"static/stress_mt{suffix}", ZNE=ZNE, dim2=True))


AVGRERR = np.array(AVGRERR)
print(AVGRERR)
print(np.mean(AVGRERR), np.min(AVGRERR), np.max(AVGRERR))
if np.mean(AVGRERR) > 0.001:
    raise ValueError

