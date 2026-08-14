# ---------------------------------------------------------------------------------
# BEGIN GRN
import numpy as np
import pygrt

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

# norths/easts 各为三个元素: start/stop/step (km)
norths = [-3.0, 3.0, 0.15]
easts = [-2.5, 2.5, 0.15]
# 可以设置 distarr 来指定震中距序列
# pymod.compute_static_grn(
#     depsrc=2.0, deprcv=0.0, distarr=np.arange(0, 10+1e-8, 0.1),
# )
# 也可以设置 norths 和 easts 来指定 north/east 网格
pymod.compute_static_grn(
    depsrc=2.0,
    deprcv=0.0,
    norths=norths,
    easts=easts,
)
static_grn = pygrt.utils.read_static_nc("stgrn.nc")
print(static_grn.keys())
# dict_keys(['dimensions', 'variables', 'attributes'])
print(list(static_grn["variables"].keys()))
# ['model', 'depsrc', 'deprcv', 'north', 'east', 'src_va', ..., 'EXZ', 'EXR', ...]
# END GRN
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN plot func
import matplotlib.pyplot as plt
from typing import Union

def plot_static(static_syn:dict, out:Union[str,None]=None):
    vars_ = static_syn["variables"]
    north = vars_["north"]["data"]
    east = vars_["east"]["data"]
    Z = vars_["Z"]["data"]
    N = vars_["N"]["data"]
    E = vars_["E"]["data"]

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    # 设计对称色标
    m = np.max(np.abs(Z)) * 1.2
    pcm = ax.pcolormesh(east, north, Z, cmap='bwr', vmin=-m, vmax=m)
    ax.quiver(east, north, E, N, angles='uv', pivot='mid')
    ax.set_ylim([north[0], north[-1]])
    ax.set_xlim([east[0], east[-1]])
    ax.set_aspect('equal')
    cbar = fig.colorbar(pcm, ax=ax, label='Z(cm)')
    cbar.formatter.set_powerlimits((0, 0))

    if out is not None:
        fig.savefig(out, bbox_inches='tight')
# END plot func
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN REUSE STGRN
# 若仅使用已算好的静态格林函数做合成，构造时只需指定 stgrn，无需 modelpath
pymod = pygrt.PyModel1D(stgrn="stgrn.nc")
# END REUSE STGRN
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN EX
static_syn = pymod.compute_static_syn(
    scale=1e24,
    output_path="stsyn_ex.nc",
    source="EX",
    zne=True,
    return_result=True,
)
print(list(static_syn["variables"].keys()))
# ['north', 'east', 'Z', 'N', 'E']
plot_static(static_syn, "syn_ex.svg")
# END SYN EX
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN SF
static_syn = pymod.compute_static_syn(
    scale=1e16,
    output_path="stsyn_sf.nc",
    source="SF",
    force=(1, -0.5, 2),
    zne=True,
    return_result=True,
)
print(list(static_syn["variables"].keys()))
# ['north', 'east', 'Z', 'N', 'E']
plot_static(static_syn, "syn_sf.svg")
# END SYN SF
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN DC
static_syn = pymod.compute_static_syn(
    scale=1e24,
    output_path="stsyn_dc.nc",
    source="DC",
    strike=33,
    dip=50,
    rake=120,
    zne=True,
    return_result=True,
)
print(list(static_syn["variables"].keys()))
# ['north', 'east', 'Z', 'N', 'E']
plot_static(static_syn, "syn_dc.svg")
# END SYN DC
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN SYN DC2
static_syn = pymod.compute_static_syn(
    scale=1e24,
    output_path="stsyn_dc2.nc",
    source="DC",
    strike=33,
    dip=90,
    rake=0,
    zne=True,
    return_result=True,
)
print(list(static_syn["variables"].keys()))
# ['north', 'east', 'Z', 'N', 'E']
plot_static(static_syn, "syn_dc2.svg")
# END SYN DC2
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN SYN TS
static_syn = pymod.compute_static_syn(
    scale=1e24,
    output_path="stsyn_ts.nc",
    source="TS",
    strike=33,
    dip=50,
    zne=True,
    return_result=True,
)
print(list(static_syn["variables"].keys()))
# ['north', 'east', 'Z', 'N', 'E']
plot_static(static_syn, "syn_ts.svg")
# END SYN TS
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN TS2
static_syn = pymod.compute_static_syn(
    scale=1e24,
    output_path="stsyn_ts2.nc",
    source="TS",
    strike=33,
    dip=90,
    zne=True,
    return_result=True,
)
print(list(static_syn["variables"].keys()))
# ['north', 'east', 'Z', 'N', 'E']
plot_static(static_syn, "syn_ts2.svg")
# END SYN TS2
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN SYN MT
static_syn = pymod.compute_static_syn(
    scale=1e24,
    output_path="stsyn_mt.nc",
    source="MT",
    moment_tensor=(0.1, -0.2, 1.0, 0.3, -0.5, -2.0),
    zne=True,
    return_result=True,
)
print(list(static_syn["variables"].keys()))
# ['north', 'east', 'Z', 'N', 'E']
plot_static(static_syn, "syn_mt.svg")
# END SYN MT
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN SYN MT2
static_syn = pymod.compute_static_syn(
    scale=1e24,
    output_path="stsyn_mt2.nc",
    source="MT",
    moment_tensor=(0, -0.2, 0, 0, 0, 0),
    zne=True,
    return_result=True,
)
print(list(static_syn["variables"].keys()))
# ['north', 'east', 'Z', 'N', 'E']
plot_static(static_syn, "syn_mt2.svg")
# END SYN MT2
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN NEW XY
static_syn = pymod.compute_static_syn(
    scale=1e24,
    output_path="stsynXY_dc2.nc",
    source="DC",
    strike=33,
    dip=90,
    rake=0,
    zne=True,
    norths=[-3.0, 3.0, 0.2],
    easts=[-2.5, 2.5, 0.25],
    return_result=True,
)
print(list(static_syn["variables"].keys()))
# ['north', 'east', 'Z', 'N', 'E']
plot_static(static_syn, "synXY_dc2.svg")
# END NEW XY
# ---------------------------------------------------------------------------------

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in [
    "stgrn.nc",
    "stsyn_ex.nc", "stsyn_sf.nc", "stsyn_dc.nc", "stsyn_dc2.nc",
    "stsyn_ts.nc", "stsyn_ts2.nc", "stsyn_mt.nc", "stsyn_mt2.nc",
    "stsynXY_dc2.nc",
]:
    p = Path(name)
    if p.is_file():
        p.unlink(missing_ok=True)
