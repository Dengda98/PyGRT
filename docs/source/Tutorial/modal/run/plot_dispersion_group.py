import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf_file

def plot_dispersion(path, title, ax):
    with netcdf_file(path, "r", mmap=False) as f:
        freq = f.variables["freq"][:]
        for i in range(f.dimensions['mode']):
            c = np.array(f.variables['u'][:, i])
            c[c == 0.0] = np.nan
            ax.plot(freq, c, lw=0.6, marker='o', ms=0, c='b', label='C')

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Group velocity (km/s)")
    ax.set_title(title)
    ax.set_xmargin(0)
    ax.set_xlim(xmin=0)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
plot_dispersion("group_R.nc", "Rayleigh", ax1)
plot_dispersion("group_L.nc", "Love", ax2)

fig.savefig("dispersion_py_group.svg")
