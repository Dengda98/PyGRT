import numpy as np
import pygrt 
import matplotlib.pyplot as plt
import glob, os

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.axes import Axes


# 定义灰色映射函数
def gray_mapping(x):
    """根据输入x计算对应的灰度值"""
    return (0.7 - x / 0.7) ** 0.8

def _plot_one(ax:Axes, pattern:str, ktype:str):
    norm = 0.0
    yLst = []
    evdpLst = []
    maxnorm = 0
    for path in glob.glob(pattern):
        data = pygrt.utils.read_statsfile(path)

        evdp = float(path.split("/")[-2].split("_")[1])

        karr = data['k']
        Farr = data[ktype]

        norm = np.linalg.norm(np.real(Farr), np.inf)
        yLst.append(np.real(Farr))
        evdpLst.append(evdp)
        maxnorm = max(maxnorm, norm)

    ax.axhline(y=0, color='k', ls='--', linewidth=1)
    
    for y, evdp in zip(yLst, evdpLst):
        ax.plot(karr, y/maxnorm * 1e8, c=str(gray_mapping(evdp)), zorder=int(evdp*1e6))
    
    ax.set_xlim(0, karr[-1])
    ax.grid(linewidth=0.4)


srctypes = ['EX', 'VF', 'HF', 'DD', 'DS', 'SS']
modrLst = [0, 0, 1, 0, 1, 2]
qwvLst = ['q', 'w', 'v']
fig = plt.figure(figsize=(14, 2.*len(srctypes)))
gs = fig.add_gridspec(len(srctypes), 3, hspace=0.4, wspace=0.2)
for isrc, src in enumerate(srctypes):
    for iqwv, qwv in enumerate(qwvLst):
        if modrLst[isrc] == 0 and qwv == 'v':
            continue
        key = f"{src}_{qwv}"

        ax = fig.add_subplot(gs[isrc, iqwv])

        _plot_one(ax, "GRN_grtstats/mod_*_0/K_0025_*", key)
        ax.set_title(key)

        ax.set_yscale('symlog')
        ax.set_yticks([-1e9, -1e6, -1e3, 0, 1e3, 1e6, 1e9])

        if isrc == len(srctypes)-1:
            ax.set_xlabel(r"$k / \text{km}^{-1}$", fontsize=12)
        else:
            ax.set_xticklabels([])

# 生成自定义灰度色图
x_vals = np.linspace(0, 0.4, 256)  # 和深度范围保持一致
gray_values = gray_mapping(x_vals)
# 创建自定义 colormap，使用灰色系
cmap = LinearSegmentedColormap.from_list(
    'custom_gray', 
    [(g, g, g) for g in gray_values]
)

# 添加色标棒
cax = fig.add_axes([0.35, 0.91, 0.37, 0.01])
sc = ax.scatter([], [], c=[], cmap=cmap, s=50)
cbar = fig.colorbar(sc, cax=cax, orientation='horizontal', 
                    ticks=[], 
                    pad=0.04, aspect=30)
cbar.ax.text(-0.02, 0.5, 'Shallow', transform=cbar.ax.transAxes, ha='right', va='center', fontsize=14)
cbar.ax.text(1.02, 0.5, 'Deep', transform=cbar.ax.transAxes, ha='left', va='center', fontsize=14)

fig.savefig("deep_shallow_kernel.svg", bbox_inches='tight')