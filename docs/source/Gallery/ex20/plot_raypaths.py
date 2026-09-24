import math

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import brentq


NU = 0.25
R = 20.0
SOURCE_DEPTH = 5.0
RECEIVER_DEPTH = 4.0
K = math.sqrt(0.5 * (1.0 - 2.0 * NU) / (1.0 - NU))
THETA_C = math.asin(K)


def conversion_split(p_depth, s_depth):
    """Return the surface conversion point for a P-S path."""

    def residual(horizontal_split):
        p_path = math.hypot(horizontal_split, p_depth)
        s_path = math.hypot(R - horizontal_split, s_depth)
        return (R - horizontal_split) / s_path - K * horizontal_split / p_path

    left = 1e-12 * (1.0 + R)
    right = R - left
    return brentq(residual, left, right)


def add_path(axis, points, phase, style, label_point, label_offset):
    x_values, depth_values = zip(*points)
    axis.plot(x_values, depth_values, color=style[0], linestyle=style[1], linewidth=2.0)
    axis.annotate(
        phase,
        xy=label_point,
        xytext=label_offset,
        textcoords="offset points",
        color=style[0],
        fontsize=11,
        ha="center",
        va="center",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5),
    )


source = (0.0, SOURCE_DEPTH)
receiver = (R, RECEIVER_DEPTH)
reflection_x = R * SOURCE_DEPTH / (SOURCE_DEPTH + RECEIVER_DEPTH)
ps_x = conversion_split(SOURCE_DEPTH, RECEIVER_DEPTH)
sp_x = conversion_split(RECEIVER_DEPTH, SOURCE_DEPTH)
sps_x1 = SOURCE_DEPTH * math.tan(THETA_C)
sps_x2 = R - RECEIVER_DEPTH * math.tan(THETA_C)

styles = {
    "P": ("#1f77b4", "-"),
    "PP": ("#2ca02c", "-"),
    "PS": ("#9467bd", "-"),
    "SP": ("#8c564b", "--"),
    "sPs": ("#e377c2", ":"),
}

fig, axis = plt.subplots(figsize=(9, 5))
axis.axhline(0.0, color="black", linewidth=3.0, zorder=1)

axis.scatter(*source, color="#d62728", s=70, zorder=5, label="Source")
axis.scatter(*receiver, color="black", s=55, zorder=5, label="Receiver")
axis.text(
    source[0], source[1] + 0.5, "Source", color="#d62728", ha="center", va="top",
    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5),
)
axis.text(
    receiver[0], receiver[1] + 0.5, "Receiver", color="black", ha="center", va="top",
    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5),
)

add_path(axis, [source, receiver], "P(S)", styles["P"], (R * 0.48, 4.55), (0, -13))

reflection_point = (reflection_x, 0.0)
add_path(
    axis,
    [source, reflection_point, receiver],
    "PP(SS)",
    styles["PP"],
    (0.5 * reflection_x, SOURCE_DEPTH * 0.28),
    (0, -12),
)

ps_point = (ps_x, 0.0)
sp_point = (sp_x, 0.0)
add_path(
    axis,
    [source, ps_point, receiver],
    "PS",
    styles["PS"],
    (0.52 * (ps_x + R), 0.5 * RECEIVER_DEPTH),
    (0, -12),
)
add_path(
    axis,
    [source, sp_point, receiver],
    "SP",
    styles["SP"],
    (0.5 * sp_x, 0.5 * SOURCE_DEPTH),
    (0, 12),
)

sps_points = [source, (sps_x1, 0.0), (sps_x2, 0.0), receiver]
add_path(axis, sps_points, "sPs", styles["sPs"], (sps_x1-1, 0.0), (0, -15))

legend_lines = [
    Line2D([0], [0], color=styles[name][0], linestyle=styles[name][1], linewidth=2.0, label=label)
    for name, label in (("P", "P(S)"), ("PP", "PP(SS)"), ("PS", "PS"), ("SP", "SP"), ("sPs", "sPs"))
]
axis.legend(handles=legend_lines, ncol=3, loc="upper right", frameon=True, framealpha=0.85, fontsize=9)
axis.set_xlim(-1.5, R + 1.5)
axis.set_ylim(max(SOURCE_DEPTH, RECEIVER_DEPTH) + 1.0, -2.0)
axis.set_xlabel("Horizontal distance (km)")
axis.set_ylabel("Depth (km)")
axis.set_title("Lamb3 phase ray paths", pad=8)
axis.grid(axis="x", color="0.88", linewidth=0.7, zorder=0)
axis.set_aspect("equal")
axis.set_axisbelow(True)
fig.tight_layout()
fig.savefig("ray_paths.svg", bbox_inches="tight")
plt.close(fig)
