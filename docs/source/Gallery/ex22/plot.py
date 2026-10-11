"""绘制《地震学中的 Lamb 问题（下）》图 5.8.3、5.8.4 的有限断层位移"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from obspy import read

plt.rcParams.update({"font.family": "DejaVu Serif", "mathtext.fontset": "stix", "font.size": 11})

for mechanism, number in [("strike_slip", "5_8_3"), ("thrust", "5_8_4")]:
    fig, axes = plt.subplots(1, 3, figsize=(10, 7), sharex=True, sharey=True, layout='constrained')
    spacing = 8.5 if mechanism == "strike_slip" else 10.0
    limits = (-15, 60) if mechanism == "strike_slip" else (-18, 74)
    scale_top = 54 if mechanism == "strike_slip" else 68

    for index, dx in enumerate(("10", "5", "2.5", "2", "1", "0.5")):
        offset = (5-index)*spacing
        # 书中 x1、x2、x3 对应 North、East、Down，因此垂直分量取 -Z
        for axis, component, sign in zip(axes, "NEZ", (1, 1, -1)):
            if dx == "0.5":
                surface = read(f"results/{mechanism}_dx{dx}_z0/{component}.sac")[0]
                axis.plot(surface.times()+surface.stats.sac.b, sign*surface.data+offset, color="0.68", lw=3.2)
            trace = read(f"results/{mechanism}_dx{dx}_z2/{component}.sac")[0]
            axis.plot(trace.times()+trace.stats.sac.b, sign*trace.data+offset, color="0.08", lw=1.5)
            axis.text(0.28, offset+1.1, rf"$\Delta x={dx}\ \mathrm{{km}}$", fontsize=11)

    for component, axis in enumerate(axes, start=1):
        axis.set_xlim(0, 12)
        axis.set_ylim(*limits)
        axis.set_xticks((0, 3, 6, 9, 12))
        axis.set_yticks([i*spacing for i in range(6)])
        axis.tick_params(direction="in", top=True, right=True, labelleft=False)
        axis.set_xlabel(r"$t\,/\,\mathrm{s}$", fontsize=13)
        axis.text(0.13, 0.04, rf"$u_{component}$", transform=axis.transAxes, fontsize=14)
        axis.plot([10.0, 9.45, 9.45, 10.0], [scale_top, scale_top, scale_top-5, scale_top-5], color="black", lw=1.0)
        axis.text(9.725, scale_top+1, "5 cm", ha="center", va="bottom", fontsize=12)

    fig.savefig(f"fig_{number}.svg")
    plt.close(fig)
