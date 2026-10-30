import matplotlib.pyplot as plt
import numpy as np
import pygrt


static_syn = pygrt.io.read_nc_variables("stsyn_points.nc")
fig, ax = pygrt.plot.plot_static_displacement(static_syn, title="Static displacement for -M33/90/0")
ax.plot(0.0, 0.0, marker="*", markersize=12, lw=0, color="black", label="Epicenter")
limit = 1.15 * max(np.max(np.abs(static_syn["north"])), np.max(np.abs(static_syn["east"])))
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.legend(loc="upper right")
fig.savefig("syn_points.svg", dpi=200, bbox_inches="tight")
plt.close(fig)
