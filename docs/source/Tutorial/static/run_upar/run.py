import matplotlib.pyplot as plt
import numpy as np
import pygrt

def plot6(data:dict, title:str, out:str|None=None):
    vars_ = data["variables"]
    norths = vars_["north"]["data"]
    easts = vars_["east"]["data"]
    chs = sorted(
        [k for k in vars_ if k.startswith(title.lower() + "_")],
        reverse=True,
    )
    fig, axs = plt.subplots(len(chs)//3, 3, figsize=(10, len(chs)))
    axs = axs.ravel()

    MAX = 0
    for i in range(len(chs)):
        ch = chs[i]
        m = np.max(np.abs(vars_[ch]["data"]))
        if m > MAX:
            MAX = m

    for i in range(len(chs)):
        ax = axs[i]
        ch = chs[i]
        arr = vars_[ch]["data"]
        vmin = vmax = None
        if np.max(np.abs(arr))/MAX < 1e-5:
            vmin = -1
            vmax = 1

        pcm = ax.pcolormesh(easts, norths, arr, shading='nearest', vmin=vmin, vmax=vmax, rasterized=True)
        ax.set_aspect('equal')
        ax.set_title(ch)
        cbar = fig.colorbar(pcm, ax=ax)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_normal(pcm)

    fig.suptitle(title)

    if out is not None:
        fig.savefig(out, bbox_inches='tight')


pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

# norths/easts 各为三个元素: start/stop/step (km)
# 传入 calc_upar=True 可计算空间导数
pymod.compute_static_grn(
    depsrc=2.0,
    deprcv=0.0,
    norths=[-3.0, 3.0, 0.15],
    easts=[-2.5, 2.5, 0.15],
    calc_upar=True,
)

# 传入 calc_upar=True 可计算空间导数
# 传入 zne=True 返回 ZNE 分量
pymod.compute_static_syn(
    scale=1e24,
    output_path="stsyn_dc_zne.nc",
    strike=33,
    dip=50,
    rake=120,
    zne=True,
    calc_upar=True,
)

# 计算应变 / 旋转 / 应力，结果写回同一 nc 文件
static_strain = pygrt.utils.compute_strain("stsyn_dc_zne.nc", return_result=True)
static_rotation = pygrt.utils.compute_rotation("stsyn_dc_zne.nc", return_result=True)
static_stress = pygrt.utils.compute_stress("stsyn_dc_zne.nc", return_result=True)

plot6(static_strain, "Strain", 'static_strain.svg')
plot6(static_rotation, "Rotation", 'static_rotation.svg')
plot6(static_stress, "Stress", 'static_stress.svg')

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in ["stgrn.nc", "stsyn_dc_zne.nc"]:
    p = Path(name)
    if p.is_file():
        p.unlink(missing_ok=True)
