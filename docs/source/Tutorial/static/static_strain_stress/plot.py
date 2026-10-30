import matplotlib.pyplot as plt
import numpy as np
import pygrt

def plot6(data: dict, title: str, out: str | None = None):
    norths = data["north"]
    easts = data["east"]
    chs = sorted([k for k in data if k.startswith(title.lower() + "_")], reverse=True)
    fig, axs = plt.subplots(len(chs)//3, 3, figsize=(10, len(chs)))
    axs = axs.ravel()

    MAX = 0
    for i in range(len(chs)):
        ch = chs[i]
        m = np.max(np.abs(data[ch]))
        if m > MAX:
            MAX = m

    for i in range(len(chs)):
        ax = axs[i]
        ch = chs[i]
        arr = data[ch]
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
    plt.close(fig)


static_tensor = pygrt.io.read_nc_variables("stsyn_dc_zne.nc")
plot6(static_tensor, "Strain", 'static_strain.svg')
plot6(static_tensor, "Rotation", 'static_rotation.svg')
plot6(static_tensor, "Stress", 'static_stress.svg')

# BEGIN COULOMB PLOT
fig, axes = pygrt.plot.plot_static_coulomb(static_tensor, outpath="static_coulomb.svg")
plt.close(fig)
# END COULOMB PLOT
