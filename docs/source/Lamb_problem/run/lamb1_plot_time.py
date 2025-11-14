# BEGIN LAMB1
import pygrt
import numpy as np

ts = np.arange(0, 2, 1e-4)
u = pygrt.utils.solve_lamb1(0.25, ts, 30)
# END LAMB1

import matplotlib.pyplot as plt

fig, axs = plt.subplots(3, 3, figsize=(10, 5), sharex=True)
for i in range(3):
    for j in range(3):
        ax = axs[i, j]
        ax.plot(ts, u[:, i, j])
        ax.set_xlim(0, 2)
        ax.set_ylim(-2, 2)

        ax.text(0.1, 0.9, rf"$\bar{{G}}^H_{{{i+1}{j+1}}}$", transform=ax.transAxes, ha='left', va='top', fontsize=12)

fig.savefig("lamb1_time.png", dpi=100, bbox_inches='tight')