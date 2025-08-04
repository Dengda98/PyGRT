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
        if np.isnan(rerr) or np.isinf(rerr):
            rerr = 0.0
        error += rerr
        nerr += 1 

    return error/nerr


def static_compare3(resDct:dict, c_prefix:str):
    """return average relative error"""

    nx = len(resDct['_xarr'])
    ny = len(resDct['_yarr'])

    def _get_c_resDct(c_prefix:str):
        with open(c_prefix, 'r') as f:
            lines = [line.strip() for line in f if line.startswith('#')]
        chs = lines[-1].split()[1:]

        c_res = np.loadtxt(c_prefix)
        c_resDct = {}
        for i in range(2, c_res.shape[1]):
            c_resDct[chs[i]] = c_res[:,i].reshape((nx, ny), order='F')

        return c_resDct
    
    c_resDct = _get_c_resDct(c_prefix)

    error = 0.0
    nerr = 0

    for k in c_resDct.keys():
        val1 = resDct[k]
        val2 = c_resDct[k]
        rerr = np.mean(np.abs(val1 - val2) / np.abs(val2))
        if np.isnan(rerr) or np.isinf(rerr):
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
    
    st = pygrt.utils.gen_syn_from_gf_EX(st_grn, S, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_triangle_wave(0.4, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, "syn_ex/cout", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    strotation = pygrt.utils.compute_rotation(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, "syn_ex/cout.strain.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(strotation, "syn_ex/cout.rotation.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, "syn_ex/cout.stress.", ZNE=ZNE, dim2=True))
    

    st = pygrt.utils.gen_syn_from_gf_SF(st_grn, S, fn, fe, fz, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_trap_wave(0.1, 0.3, 0.6, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, "syn_sf/cout", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    strotation = pygrt.utils.compute_rotation(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, "syn_sf/cout.strain.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(strotation, "syn_sf/cout.rotation.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, "syn_sf/cout.stress.", ZNE=ZNE, dim2=True))

    st = pygrt.utils.gen_syn_from_gf_DC(st_grn, S, stk, dip, rak, az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_parabola_wave(0.6, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, "syn_dc/cout", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    strotation = pygrt.utils.compute_rotation(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, "syn_dc/cout.strain.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(strotation, "syn_dc/cout.rotation.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, "syn_dc/cout.stress.", ZNE=ZNE, dim2=True))

    st = pygrt.utils.gen_syn_from_gf_MT(st_grn, S, [M11,M12,M13,M22,M23,M33], az, ZNE=ZNE, calc_upar=True)
    sigs = pygrt.sigs.gen_ricker_wave(3, dt)
    pygrt.utils.stream_convolve(st, sigs)
    AVGRERR.append(compare3(st, "syn_mt/cout", ZNE=ZNE))
    ststrain = pygrt.utils.compute_strain(st)
    strotation = pygrt.utils.compute_rotation(st)
    ststress = pygrt.utils.compute_stress(st)
    AVGRERR.append(compare3(ststrain, "syn_mt/cout.strain.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(strotation, "syn_mt/cout.rotation.", ZNE=ZNE, dim2=True))
    AVGRERR.append(compare3(ststress, "syn_mt/cout.stress.", ZNE=ZNE, dim2=True))


#-------------------------- Static -----------------------------------------
xarr = np.arange(-3, 3.1, 0.6)
yarr = np.arange(-4, 4.1, 0.8)
static_grn = pymod.compute_static_grn(xarr, yarr, calc_upar=True)

for ZNE in [False, True]:
    suffix = "-N" if ZNE else ""
    static_syn = pygrt.utils.gen_syn_from_gf_EX(static_grn, S, ZNE=ZNE, calc_upar=True)
    AVGRERR.append(static_compare3(static_syn, f"static/stsyn_ex{suffix}"))
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    AVGRERR.append(static_compare3(ststrain, f"static/strain_ex{suffix}"))
    AVGRERR.append(static_compare3(strotation, f"static/rotation_ex{suffix}"))
    AVGRERR.append(static_compare3(ststress, f"static/stress_ex{suffix}"))
    print(static_syn, ststrain, strotation, ststress)

    static_syn = pygrt.utils.gen_syn_from_gf_SF(static_grn, S, fn, fe, fz, ZNE=ZNE, calc_upar=True)
    AVGRERR.append(static_compare3(static_syn, f"static/stsyn_sf{suffix}"))
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    AVGRERR.append(static_compare3(ststrain, f"static/strain_sf{suffix}"))
    AVGRERR.append(static_compare3(strotation, f"static/rotation_sf{suffix}"))
    AVGRERR.append(static_compare3(ststress, f"static/stress_sf{suffix}"))
    print(static_syn, ststrain, strotation, ststress)

    static_syn = pygrt.utils.gen_syn_from_gf_DC(static_grn, S, stk, dip, rak, ZNE=ZNE, calc_upar=True)
    AVGRERR.append(static_compare3(static_syn, f"static/stsyn_dc{suffix}"))
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    AVGRERR.append(static_compare3(ststrain, f"static/strain_dc{suffix}"))
    AVGRERR.append(static_compare3(strotation, f"static/rotation_dc{suffix}"))
    AVGRERR.append(static_compare3(ststress, f"static/stress_dc{suffix}"))
    print(static_syn, ststrain, strotation, ststress)

    static_syn = pygrt.utils.gen_syn_from_gf_MT(static_grn, S, [M11,M12,M13,M22,M23,M33], ZNE=ZNE, calc_upar=True)
    AVGRERR.append(static_compare3(static_syn, f"static/stsyn_mt{suffix}"))
    ststrain = pygrt.utils.compute_strain(static_syn)
    strotation = pygrt.utils.compute_rotation(static_syn)
    ststress = pygrt.utils.compute_stress(static_syn)
    AVGRERR.append(static_compare3(ststrain, f"static/strain_mt{suffix}"))
    AVGRERR.append(static_compare3(strotation, f"static/rotation_mt{suffix}"))
    AVGRERR.append(static_compare3(ststress, f"static/stress_mt{suffix}"))
    print(static_syn, ststrain, strotation, ststress)


AVGRERR = np.array(AVGRERR)
print(AVGRERR)
print(np.mean(AVGRERR), np.min(AVGRERR), np.max(AVGRERR))
if np.mean(AVGRERR) > 0.01:
    raise ValueError

