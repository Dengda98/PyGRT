"""绘制单个地表接收点的动态位移和投影应力"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from obspy import read


plt.rcParams.update({"font.size": 11, "svg.fonttype": "none"})

for name, channels, labels, factor, unit in (
    ("displacement", ("Z", "N", "E"), ("Up", "North", "East"), 1.0, "Displacement (cm)"),
    ("stress", ("sigma_n", "tau_s", "coulomb"),
     (r"Normal stress $\Delta\sigma_n$", r"Shear stress $\Delta\tau_s$", r"Coulomb stress $\Delta CFS$"), 1e-4, "Stress (kPa)"),
):
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True, layout="constrained")
    for ax, channel, label, color in zip(axes, channels, labels, ("C0", "C1", "C3")):
        trace = read(f"syn/{channel}.sac")[0]
        time = trace.times() + trace.stats.sac.b
        ax.plot(time, trace.data * factor, color=color, lw=1.4, label=label)
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_ylabel(unit)
        ax.set_xlim(0, 400)
        ax.grid(alpha=0.25, lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(loc="upper right", fontsize=11, frameon=False)
    axes[-1].set_xlabel("Time since rupture onset (s)")
    fig.savefig(f"{name}.svg")
    plt.close(fig)
