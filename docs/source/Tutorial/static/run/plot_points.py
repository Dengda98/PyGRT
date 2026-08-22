import numpy as np
import pygrt
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def plot_static_points(static_syn, output_path):
    variables = static_syn["variables"]
    north = variables["north"]["data"]
    east = variables["east"]["data"]
    z = variables["Z"]["data"]
    n = variables["N"]["data"]
    e = variables["E"]["data"]

    zmax = max(float(np.max(np.abs(z))), np.finfo(float).eps)

    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
    points = ax.scatter(
        east, north,
        s=180,
        c=z,
        cmap="turbo",
        norm=Normalize(-zmax, zmax),
        linewidths=0.,
        zorder=10,
    )
    ax.quiver(
        east, north, e, n,
        color="black", angles="xy", scale_units="xy", scale=None, width=0.003,
        zorder=20,
    )
    ax.plot(0.0, 0.0, marker="*", markersize=12, lw=0, color="black", label="Epicenter")

    limit = 1.15 * max(np.max(np.abs(north)), np.max(np.abs(east)))
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal")
    ax.set_xlabel("East (km)")
    ax.set_ylabel("North (km)")
    ax.set_title("Static displacement for -M33/90/0")
    ax.legend(loc="upper right")
    ax.grid(lw=0.4)
    fig.colorbar(points, ax=ax, label="Vertical displacement Z (cm)", shrink=0.7)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


static_syn = pygrt.utils.read_static_nc("stsyn_points.nc")
plot_static_points(static_syn, "syn_points.svg")