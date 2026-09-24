import matplotlib.pyplot as plt

# BEGIN LAMB1
import pygrt
import numpy as np

tbar = np.arange(0, 2, 1e-4)
u = pygrt.utils.lamb1(nu=0.25, tbar=tbar, azimuth=30)
# END LAMB1

fig, axs = plt.subplots(3, 3, figsize=(10, 5), sharex=True)
for i in range(3):
    for j in range(3):
        ax = axs[i, j]
        ax.plot(tbar, u[:, i, j])
        ax.set_xlim(0, 2)
        ax.set_ylim(-2, 2)

        ax.text(0.1, 0.9, rf"$\bar{{G}}^H_{{{i+1}{j+1}}}$", transform=ax.transAxes, ha='left', va='top', fontsize=12)

fig.savefig("lamb1_time.svg", bbox_inches='tight')
plt.close()

# BEGIN LAMB1 MOVING
import pygrt
import numpy as np

tbar = np.arange(0, 2, 1e-4)
u = pygrt.utils.lamb1(nu=0.25, tbar=tbar, azimuth=30, cbar=0.1)
# END LAMB1 MOVING

fig, axs = plt.subplots(1, 3, figsize=(10, 2))
for i in range(3):
    ax = axs[i]
    ax.plot(tbar, u[:, i])
    ax.set_xlim(0, 2)
    ax.set_ylim(-2, 2)

    ax.text(0.1, 0.9, rf"$\bar{{u}}{{{i+1}}}$", transform=ax.transAxes, ha='left', va='top', fontsize=12)

fig.savefig("lamb1_time_moving.svg", bbox_inches='tight')