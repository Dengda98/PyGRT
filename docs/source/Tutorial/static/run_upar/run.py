import numpy as np
import pygrt 

modarr = np.loadtxt("milrow")

pymod = pygrt.PyModel1D(modarr, depsrc=2.0, deprcv=0.0)

xarr = np.linspace(-3, 3, 41)
yarr = np.linspace(-2.5, 2.5, 33)
# 传入calc_upar=True可计算空间导数
static_grn = pymod.compute_static_grn(xarr, yarr, calc_upar=True)

# 传入calc_upar=True可计算空间导数
# 传入ZNE=True返回ZNE分量
static_syn = pygrt.utils.gen_syn_from_gf_DC(static_grn, M0=1e24, strike=33, dip=50, rake=120, ZNE=True, calc_upar=True)

# 计算应变
static_strain = pygrt.utils.compute_strain(static_syn)

# 计算应力
static_stress = pygrt.utils.compute_stress(static_syn)


# 简单绘制结果
import matplotlib.pyplot as plt

def plot6(data:dict, title:str, out:str|None=None):
    chs = ['ZZ', 'ZN', 'ZE', 'NN', 'NE', 'EE']
    xarr = data['_xarr']
    yarr = data['_yarr']
    fig, axs = plt.subplots(2, 3, figsize=(10, 6))
    axs = axs.ravel()

    MAX = 0
    for i in range(6):
        ch = chs[i]
        m = np.max(np.abs(data[ch]))
        if m > MAX:
            MAX = m

    for i in range(6):
        ax = axs[i]
        ch = chs[i]
        vmin = vmax = None
        if np.max(np.abs(data[ch]))/MAX < 1e-5:
            vmin = -1
            vmax = 1

        pcm = ax.pcolormesh(yarr, xarr, data[ch], shading='nearest', vmin=vmin, vmax=vmax)
        ax.set_title(ch)
        cbar = fig.colorbar(pcm, ax=ax)
        cbar.formatter.set_powerlimits((0, 0))
        cbar.update_normal(pcm)

    fig.suptitle(title)

    if out is not None:
        fig.tight_layout()
        fig.savefig(out, dpi=100)

plot6(static_strain, "Strain", 'static_strain.png')
plot6(static_stress, "Stress", 'static_stress.png')