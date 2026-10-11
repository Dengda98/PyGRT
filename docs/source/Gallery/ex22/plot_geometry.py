"""调用断层绘图函数，展示两种震源断层及接收点的位置"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pygrt

for mechanism, title in (("strike_slip", "Left-lateral strike-slip fault"), ("thrust", "Reverse fault")):
    fig, (ax_map, ax_3d) = pygrt.utils.plot_sCoulomb_faults_geometry(f"{mechanism}.inr")
    fig.set_size_inches(12, 6)
    figure_title = fig.suptitle(title, fontsize=14)

    # 接收点的水平位置重合，sCoulomb 绘图坐标为 East、North、Depth
    ax_map.plot(0, 20, "kv", ms=7, label="Receivers (depth 0 / 2 km)", zorder=5)
    ax_3d.scatter(0, 20, 2, color="black", marker="v", s=50, depthshade=False, label="Receiver (depth 2 km)")
    ax_3d.scatter(0, 20, 0, color="0.55", marker="o", s=40, depthshade=False, label="Receiver (surface)")

    # 两张图使用相同的坐标范围，便于比较断层与接收点的位置
    for axis in (ax_map, ax_3d):
        axis.set_xlim(-13, 13)
        axis.set_ylim(-12, 24)
    ax_3d.set_zlim(11.2, 0)
    ax_3d.set_zticks((0, 2, 4, 6, 8, 10))
    ax_3d.set_aspect("equal")
    ax_map.legend(loc="lower left", fontsize=8)
    ax_3d.legend(loc="upper left", fontsize=8)

    fig.savefig(f"geometry_{mechanism}.svg", bbox_inches="tight", bbox_extra_artists=(figure_title, ax_3d.zaxis.label))
    plt.close(fig)
