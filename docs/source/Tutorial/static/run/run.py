# BEGIN GRN
import numpy as np
import pygrt 

modarr = np.loadtxt("milrow")

pymod = pygrt.PyModel1D(modarr, depsrc=2.0, deprcv=0.0)

xarr = np.linspace(-3, 3, 41)
yarr = np.linspace(-2.5, 2.5, 33)
static_grn = pymod.compute_static_grn(xarr, yarr)
print(static_grn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'EXZ', 'VFZ', 'DDZ', 'HFZ', 'DSZ', 'SSZ', 'EXR', 'VFR', 'DDR', 'HFR', 'DSR', 'SSR', 'HFT', 'DST', 'SST'])
# END GRN

# BEGIN SYN EXP
static_syn = pygrt.utils.gen_syn_from_gf_EXP(static_grn, M0=1e24, ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'EXZ', 'EXN', 'EXE'])
# END SYN EXP


# BEGIN SYN SF
static_syn = pygrt.utils.gen_syn_from_gf_SF(static_grn, S=1e16, fN=1, fE=-0.5, fZ=2, ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'SFZ', 'SFN', 'SFE'])
# END SYN SF

# BEGIN SYN DC
static_syn = pygrt.utils.gen_syn_from_gf_DC(static_grn, M0=1e24, strike=33, dip=50, rake=120, ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'DCZ', 'DCN', 'DCE'])
# END SYN DC

# BEGIN SYN MT
static_syn = pygrt.utils.gen_syn_from_gf_MT(static_grn, M0=1e24, MT=[0.1,-0.2,1.0,0.3,-0.5,-2.0], ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'MTZ', 'MTN', 'MTE'])
# END SYN MT
