import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file
from matplotlib.lines import Line2D

def plot_eigenfunction(path, title, scale, ax):
    solidkwargs = dict(c='k', ls='-', lw=0.5)
    dashkwargs = dict(c='k', ls='--', lw=0.5)
    isRayl = False
    with netcdf_file(path, "r", mmap=False) as f:
        modes = f.variables['mode'][:]
        zarr = f.variables['z'][:]

        for i in range(f.dimensions['mode']):
            offset = modes[i]
            data = np.array(f.variables['eigfn'][0, i, :, :])
            
            # 垂向索引 和 径向索引
            idx1, idx2 = 2, 0

            # Rayleigh
            if data.shape[1] >= 8:
                trace = data[:, idx1]
                norm = np.max(np.abs(trace))
                trace /= norm * scale
                ax.plot(offset + trace, zarr, **solidkwargs)
                trace = data[:, idx2]
                trace /= norm
                ax.plot(offset + trace, zarr, **dashkwargs)
                isRayl = True
            # Love
            else:
                trace = np.real(data[:, idx2])
                norm = np.max(np.abs(trace))
                trace /= norm * scale
                ax.plot(offset + trace, zarr, **solidkwargs)

        ax.yaxis.set_inverted(True)
        ax.grid()
        ax.set_ymargin(0)
        ax.set_ylabel("Depth (km)")
        ax.set_xlabel("Order of modes")
        ax.set_title(title)

        if isRayl:
            ax.legend([Line2D([],[], **solidkwargs), Line2D([],[], **dashkwargs)], 
                      ['Vertical', 'Radial'], loc='lower right', ncol=1)
        else:
            ax.legend([Line2D([],[], **solidkwargs)], 
                      ['Transverse'], loc='lower right', ncol=1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
plot_eigenfunction("egn_R.nc", "Rayleigh", 2, ax1)
plot_eigenfunction("egn_L.nc", "Love", 2, ax2)

fig.savefig("eigenfunction.svg")
