"""绘制六种断层剖分，复现《地震学中的 Lamb 问题（下）》图 5.8.2 的网格示意"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

plt.rcParams.update({"font.family": "DejaVu Serif", "mathtext.fontset": "stix", "font.size": 11})

# 两种机制使用相同的局部网格，这里读取走滑断层的实际剖分，展开到断层面
parent = np.loadtxt("strike_slip.inr", skiprows=2)
origin = parent[1:3]
along = parent[3:5] - origin
length = np.linalg.norm(along)
along /= length
sind = np.sin(np.deg2rad(parent[8]))
width = (parent[10] - parent[9]) / sind

fig, axes = plt.subplots(2, 3, figsize=(10, 4.9), layout="constrained")
for axis, dx in zip(axes.flat, (10, 5, 2.5, 2, 1, 0.5)):
    patches = np.loadtxt(f"inputs/strike_slip_dx{dx:g}.inp", skiprows=2, usecols=range(11), ndmin=2)
    starts = (patches[:, 1:3] - origin) @ along
    ends = (patches[:, 3:5] - origin) @ along
    tops = (patches[:, 9] - parent[9]) / sind
    bots = (patches[:, 10] - parent[9]) / sind

    cells = [[(start, top), (end, top), (end, bot), (start, bot)] for start, end, top, bot in zip(starts, ends, tops, bots)]
    axis.add_collection(PolyCollection(cells, facecolors="white", edgecolors="0.45", linewidths=0.5))
    axis.scatter((starts + ends) / 2, (tops + bots) / 2, s=5, color="0.15", zorder=3)
    axis.set(xlim=(0, length), ylim=(width, 0), aspect="equal", xticks=[], yticks=[])
    axis.set_title(rf"$\Delta x={dx:g}\ \mathrm{{km}}$   ({round(length/dx)} $\times$ {round(width/dx)})")

fig.supxlabel("Along strike: 20 km", fontsize=12)
fig.supylabel("Down dip: 10 km", fontsize=12)
fig.savefig("fig_5_8_2.svg")
plt.close(fig)
