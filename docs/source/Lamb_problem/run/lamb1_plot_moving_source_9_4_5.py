# BEGIN LAMB1 MOVING SOURCE 9.4.5
import matplotlib.pyplot as plt
import numpy as np
import pygrt


alpha = 1.8  # km/s
beta = 0.5  # km/s
speed_physical = 300.0 / 3600.0  # km/s
velocity_ratio = beta / alpha
nu = (1.0 - 2.0 * velocity_ratio**2) / (2.0 * (1.0 - velocity_ratio**2))
cbar = speed_physical / beta
times = np.arange(0.0, 100.0, 0.02)
display_scales = (2.0, 2.0, 2.0)


def solve_station(x1, x2):
    distance = np.hypot(x1, x2)
    azimuth = np.degrees(np.arctan2(x2, x1)) % 360.0
    tbar = beta * times / distance
    # lamb1 返回 pi^2*mu*r*u，图 9.4.5 使用 pi^2*mu*u
    return pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth, cbar=cbar) / distance


def plot_line(ax, x1, x2s, line_index, component):
    offsets = np.arange(len(x2s)-1, -1, -1, dtype=float)
    scale = display_scales[component]
    for i, (x2, offset) in enumerate(zip(x2s, offsets), start=1):
        displacement = solve_station(x1, x2)[:, component]
        ax.hlines(offset, 0.0, times[-1], color="0.65", lw=0.5, ls="-.")
        ax.plot(times, offset + scale * displacement, color="black", lw=0.8)
        ax.text(-2.0, offset, rf"$l^{{({line_index})}}_{{{i}}}$", ha="right", va="center")

        if x1 > 0.0:
            t_projection = x1 / speed_physical
            t_pole = (x1 * x1 + x2 * x2) / (x1 * speed_physical)
            ax.plot(t_projection, offset, marker="+", color="black", ms=5, mew=0.8)
            ax.plot(t_pole, offset, marker="x", color="black", ms=4, mew=0.8)


fig, axes = plt.subplots(2, 3, figsize=(9.0, 8.0), sharey=True)
x2s = np.arange(-2.5, 3.0, 1.0)
for component, name in enumerate((r"$\bar{u}_1$", r"$\bar{u}_2$", r"$\bar{u}_3$")):
    plot_line(axes[0, component], -3.0, x2s, 3, component)
    plot_line(axes[1, component], 3.0, x2s, 4, component)
    axes[0, component].text(0.08, 0.93, name, transform=axes[0, component].transAxes)
    axes[1, component].text(0.08, 0.93, name, transform=axes[1, component].transAxes)

for ax in axes.flat:
    ax.set_xlim(-15.0, 100.0)
    ax.set_ylim(-1.0, 6.0)
    ax.set_xticks((0, 25, 50, 75, 100))
    ax.set_yticks(())
    ax.set_xlabel(r"$t\,/\,\mathrm{s}$")
    ax.tick_params(direction="in", top=True, right=True)

fig.tight_layout(w_pad=1.2, h_pad=1.8)
fig.savefig("lamb1_moving_source_9_4_5.svg", bbox_inches="tight")
# END LAMB1 MOVING SOURCE 9.4.5
