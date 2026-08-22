import matplotlib.pyplot as plt
import numpy as np
import pygrt

def plot6(data:dict, title:str, out:str|None=None):
    vars_ = data["variables"]
    norths = vars_["north"]["data"]
    easts = vars_["east"]["data"]
    chs = sorted([k for k in vars_ if k.startswith(title.lower() + "_")], reverse=True)
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



# ---------------------------------------------------------------
# BEGIN
# 重新计算格林函数，因为以下计算需要位移空间导数
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

# 计算格林函数
# 传入 calc_upar=True 可计算位移空间导数
pymod.compute_static_grn(depsrc=2.0, deprcv=0.0, dists=np.arange(0.0, 10.0 + 1e-8, 0.15), calc_upar=True)

norths = [-3.0, 3.0, 0.15]
easts = [-2.5, 2.5, 0.15]

# 传入 calc_upar=True 可合成位移空间导数
# 传入 zne=True 返回 ZNE 分量
pymod.compute_static_syn(
    norths=norths, easts=easts,
    scale=1e24, strike=33, dip=50, rake=120, output_path="stsyn_dc_zne.nc",
    zne=True, calc_upar=True
)

# 计算应变 / 旋转 / 应力，结果写回同一 nc 文件
pygrt.utils.compute_strain("stsyn_dc_zne.nc")
pygrt.utils.compute_rotation("stsyn_dc_zne.nc")
pygrt.utils.compute_stress("stsyn_dc_zne.nc")
# END
# ---------------------------------------------------------------


static_tensor = pygrt.utils.read_static_nc("stsyn_dc_zne.nc")
plot6(static_tensor, "Strain", 'static_strain.svg')
plot6(static_tensor, "Rotation", 'static_rotation.svg')
plot6(static_tensor, "Stress", 'static_stress.svg')

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in ["stsyn_dc_zne.nc"]:
    p = Path(name)
    if p.is_file():
        p.unlink(missing_ok=True)
