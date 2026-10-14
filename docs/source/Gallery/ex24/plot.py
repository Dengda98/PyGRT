"""绘制四种合成方式的阶跃剪切响应，精确节点为宽灰色底层曲线"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from obspy.io.sac import SACTrace


plt.rcParams.update({"font.size": 11, "svg.fonttype": "none"})

# 精确节点先画在底层，P 波对齐结果加粗并放在最上层
methods = [
    ("Exact node",       "syn_exact",     "0.75",    "-",  4.5, 1),
    ("Linear",           "syn_linear",    "#d55e00", "-",  1.4, 2),
    ("Nearest",          "syn_nearest",   "#0072b2", "--", 1.4, 3),
    ("P-aligned linear", "syn_p_aligned", "#009e73", "-",  2.5, 4),
]
reference = SACTrace.read("syn_exact/Z.sac")
tp, ts = reference.t0, reference.t1
windows = ((0, 18), (max(0, tp - 0.5), ts + 2.5))
fig, axes = plt.subplots(3, 2, figsize=(13.5, 9.4), sharex="col", layout="constrained")

for row, component in enumerate("ZRT"):
    for label, directory, color, style, width, order in methods:
        if directory == "syn_linear":
            corners = [SACTrace.read(f"{directory}/{depth}_{distance}/{component}.sac")
                       for depth in (4.5, 5.5) for distance in (9, 11)]
            # 各角点的起点、采样间隔和点数相同，目标在网格中心，权重均为 1/4
            trace = corners[0]
            time = trace.b + np.arange(trace.npts) * trace.delta
            data = np.mean([corner.data for corner in corners], axis=0, dtype=np.float64)
        else:
            trace = SACTrace.read(f"{directory}/{component}.sac")
            time = trace.b + np.arange(trace.npts) * trace.delta
            data = trace.data
        for column, window in enumerate(windows):
            ax = axes[row, column]
            ax.plot(time, data, label=label, color=color, linestyle=style, linewidth=width, zorder=order)
            ax.set_xlim(window)

    for column in range(2):
        ax = axes[row, column]
        ax.axvline(tp, color="0.55", linewidth=0.8, zorder=0)
        ax.axvline(ts, color="0.55", linewidth=0.8, linestyle="--", zorder=0)
        ax.set_axisbelow(True)
        ax.grid(alpha=0.18)
        ax.set_ylabel(f"{component} displacement (cm)")
        ax.spines[["top", "right"]].set_visible(False)
        if row == 0:
            ax.set_title("Full waveform" if column == 0 else "Arrival window", fontsize=12)
            ax.text(tp, 0.97, " P", transform=ax.get_xaxis_transform(), va="top", color="0.35")
            ax.text(ts, 0.97, " S", transform=ax.get_xaxis_transform(), va="top", color="0.35")
        if row == 2:
            ax.set_xlabel("Time since origin (s)")

    # 红框对应右列的坐标范围，左列略增留白以完整显示边框
    xmin, xmax = axes[row, 1].get_xlim()
    ymin, ymax = axes[row, 1].get_ylim()
    axes[row, 0].margins(y=0.08)
    axes[row, 0].add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
                                   fill=False, edgecolor="red", linewidth=1.2, zorder=5))

axes[0, 0].legend(loc="upper right", framealpha=0.95, fontsize=10)
fig.savefig("comparison.svg")
plt.close(fig)
