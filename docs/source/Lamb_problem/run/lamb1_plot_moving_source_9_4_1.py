# BEGIN LAMB1 MOVING SOURCE
import matplotlib.pyplot as plt
import numpy as np
import pygrt


alpha = 1.8  # km/s
beta = 0.5  # km/s
speed_physical = 0.1e-3  # km/s
station_x1 = 3.0  # km
station_x2 = 0.5  # km

distance = np.hypot(station_x1, station_x2)
azimuth = np.degrees(np.arctan2(station_x2, station_x1))
velocity_ratio = beta / alpha
nu = (1.0 - 2.0 * velocity_ratio**2) / (2.0 * (1.0 - velocity_ratio**2))
cbar = speed_physical / beta


def solve(t):
    tbar = beta * t / distance
    fixed = pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth)[:, :, 2]
    moving = pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth, cbar=cbar)
    return fixed, moving


t_long = np.arange(0.0, 80.0, 0.08)
t_short = np.arange(0.0, 10.0, 0.008)
fixed_long, moving_long = solve(t_long)
fixed_short, moving_short = solve(t_short)

relative_error = np.linalg.norm(moving_long - fixed_long) / np.linalg.norm(fixed_long)
print(f"relative error for c/beta={cbar:.1e}: {relative_error:.3e}")
if relative_error > 1e-2:
    raise RuntimeError("The low-cbar moving-source solution does not approach the fixed-source solution")
# END LAMB1 MOVING SOURCE


fig, axes = plt.subplots(2, 3, figsize=(11, 6), sharey="row")
component_names = [r"$\bar{u}_1$", r"$\bar{u}_2$", r"$\bar{u}_3$"]

for component, name in enumerate(component_names):
    ax = axes[0, component]
    ax.plot(t_long, fixed_long[:, component], color="0.65", lw=2.2)
    ax.plot(t_long, moving_long[:, component], color="black", lw=0.8)
    ax.text(0.05, 0.88, name, transform=ax.transAxes)
    ax.set_xlim(0.0, 80.0)
    ax.set_xticks([0, 20, 40, 60, 80])

    ax = axes[1, component]
    ax.plot(t_short, fixed_short[:, component], color="0.65", lw=2.2)
    ax.plot(t_short, moving_short[:, component], color="black", lw=0.8)
    ax.text(0.05, 0.88, name, transform=ax.transAxes)
    ax.set_xlim(0.0, 10.0)
    ax.set_xticks([0, 2, 4, 6, 8, 10])

for ax in axes.flat:
    ax.set_ylim(-2.0, 2.0)
    ax.set_yticks([-2, 0, 2])
    ax.set_xlabel(r"$t\,/\,\mathrm{s}$")
    ax.tick_params(direction="in", top=True, right=True)

fig.tight_layout(w_pad=1.5, h_pad=1.5)
fig.savefig("lamb1_moving_source_9_4_1.svg", bbox_inches="tight")
