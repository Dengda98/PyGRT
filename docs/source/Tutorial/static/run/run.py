# ---------------------------------------------------------------------------------
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
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN plot func
import matplotlib.pyplot as plt
from typing import Union

def plot_static(static_syn:dict, out:Union[str,None]=None):
    fig, ax = plt.subplots(1, 1, figsize=(10,8))
    # 设计对称色标
    m = np.max(np.abs(static_syn[f'Z'])) * 1.2
    pcm = ax.pcolormesh(yarr, xarr, static_syn[f'Z'], cmap='bwr', vmin=-m, vmax=m)
    ax.quiver(yarr, xarr, static_syn[f'E'], static_syn[f'N'], 
            angles='uv', pivot='mid')
    ax.set_ylim([xarr[0], xarr[-1]])
    ax.set_xlim([yarr[0], yarr[-1]])
    ax.set_aspect('equal')
    cbar = fig.colorbar(pcm, ax=ax, label='Z(cm)')
    cbar.formatter.set_powerlimits((0, 0))

    if out is not None:
        fig.savefig(out, bbox_inches='tight')
# END plot func
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN EX
static_syn = pygrt.utils.gen_syn_from_gf_EX(static_grn, M0=1e24, ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'Z', 'N', 'E'])
plot_static(static_syn, "syn_ex.svg")
# END SYN EX
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN SF
static_syn = pygrt.utils.gen_syn_from_gf_SF(static_grn, S=1e16, fN=1, fE=-0.5, fZ=2, ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'Z', 'N', 'E'])
plot_static(static_syn, "syn_sf.svg")
# END SYN SF
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN DC
static_syn = pygrt.utils.gen_syn_from_gf_DC(static_grn, M0=1e24, strike=33, dip=50, rake=120, ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'Z', 'N', 'E'])
plot_static(static_syn, "syn_dc.svg")
# END SYN DC
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN SYN DC2
static_syn = pygrt.utils.gen_syn_from_gf_DC(static_grn, M0=1e24, strike=33, dip=90, rake=0, ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'Z', 'N', 'E'])
plot_static(static_syn, "syn_dc2.svg")
# END SYN DC2
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN SYN TS
static_syn = pygrt.utils.gen_syn_from_gf_TS(static_grn, M0=1e24, strike=33, dip=50, ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'Z', 'N', 'E'])
plot_static(static_syn, "syn_ts.svg")
# END SYN TS
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN TS2
static_syn = pygrt.utils.gen_syn_from_gf_TS(static_grn, M0=1e24, strike=33, dip=90, ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'Z', 'N', 'E'])
plot_static(static_syn, "syn_ts2.svg")
# END SYN TS2
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN MT
static_syn = pygrt.utils.gen_syn_from_gf_MT(static_grn, M0=1e24, MT=[0.1,-0.2,1.0,0.3,-0.5,-2.0], ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'Z', 'N', 'E'])
plot_static(static_syn, "syn_mt.svg")
# END SYN MT
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN SYN MT2
static_syn = pygrt.utils.gen_syn_from_gf_MT(static_grn, M0=1e24, MT=[0,-0.2,0,0,0,0], ZNE=True)
print(static_syn.keys())
# dict_keys(['_xarr', '_yarr', '_src_va', '_src_vb', '_src_rho', '_rcv_va', '_rcv_vb', '_rcv_rho', 'Z', 'N', 'E'])
plot_static(static_syn, "syn_mt2.svg")
# END SYN MT2
# ---------------------------------------------------------------------------------
