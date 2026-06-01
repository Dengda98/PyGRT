import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf_file

plt.rcParams.update({
    "mathtext.fontset": "cm"
})

def _plot(f, char, idx, ax):
    # 仅基阶
    imode = 0
    freqs = f.variables["freq"][:]
    zs = f.variables['z'][:]

    hLst = []
    for ip in range(len(freqs)):
        data = f.variables[f"{char.lower()}sens"][ip, imode, :, idx]
        freq = freqs[ip]
        h, = ax.plot(data, zs, lw=0.3, marker='o', ms=2, label=f'$f={freq:.1f} Hz$')
        hLst.append(h)

    ax.set_ymargin(0)
    ax.yaxis.set_inverted(True)
    return hLst


def plot_sensitivity(path, char):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(5,6), layout='constrained', sharey=True)

    with netcdf_file(path, "r", mmap=False) as f:
        _plot(f, char, 0, ax1)
        ax1.set_xlabel(rf"$\dfrac{{\alpha}}{{{char}}} \dfrac{{\partial {char}}}{{\partial \alpha}}$", fontsize=14)

        hLst = _plot(f, char, 1, ax2)
        ax2.set_xlabel(rf"$\dfrac{{\beta}}{{{char}}} \dfrac{{\partial {char}}}{{\partial \beta}}$", fontsize=14)
        ax2.legend(handles=hLst, loc='lower right', ncol=1)

        _plot(f, char, 2, ax3)
        ax3.set_xlabel(rf"$\dfrac{{\rho}}{{{char}}} \dfrac{{\partial {char}}}{{\partial \rho}}$", fontsize=14)

    
    return fig, (ax1, ax2, ax3)

for char in ['U', 'c']:
    fig, _ = plot_sensitivity(f"{char.lower()}sens_R.nc", char)
    fig.suptitle("Rayleigh")
    fig.savefig(f"{char.lower()}sens_R.svg")
    plt.close()

    fig, _ = plot_sensitivity(f"{char.lower()}sens_L.nc", char)
    fig.suptitle("Love")
    fig.savefig(f"{char.lower()}sens_L.svg")
    plt.close()
