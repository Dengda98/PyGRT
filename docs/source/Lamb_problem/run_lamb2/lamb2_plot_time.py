"""使用广义闭合解复现图 7.4.6、7.4.10、7.4.11 和 7.4.12"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

import numpy as np
import pygrt

"""
# BEGIN LAMB2
import numpy as np
import pygrt

tbar = np.arange(0, 2, 0.002)
epicentral_distance = 10.0
depsrc = 5.0
azimuth = 30.0
G, dG_source, dG_receiver = pygrt.utils.lamb2(
    nu=0.25, tbar=tbar, R=epicentral_distance, depsrc=depsrc, azimuth=azimuth
)
# END LAMB2
"""

EPICENTRAL_DISTANCE = 10.0
NU = 0.25
AZIMUTH = 30.0
SOURCE_DEPTHS = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
TMAX = 2.0
DTBAR = 0.002

# 图 7.4.6 中用纵向平移区分不同震源深度
GREEN_SCALE = 0.4
GREEN_OFFSETS = np.array([1.25, 0.75, 0.25, -0.25, -0.75, -1.25])
GREEN_YLIM = (-1.8, 2.2)

# 一阶空间导数的振幅较大，分别采用适合三组导数的显示缩放
DERIVATIVE_SCALES = {1: 1e-2, 2: 1e-2, 3: 1e-2}
DERIVATIVE_OFFSETS = np.array([0.25, 0.15, 0.05, -0.05, -0.15, -0.25]) * 1.4
DERIVATIVE_YLIM = (-0.7, 1.3)

COMPONENTS = [(i, j, f"{i + 1}{j + 1}") for i in range(3) for j in range(3)]


def component_label(component: str, derivative_coordinate: int | None = None) -> str:
    if derivative_coordinate is not None:
        component = f"{component},{derivative_coordinate}'"
    return rf"$\bar{{G}}^H_{{{component}}}$"


def compute_lamb2_results() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ts = np.arange(int(round(TMAX / DTBAR)) + 1, dtype=float) * DTBAR
    greens = np.empty((len(SOURCE_DEPTHS), len(ts), 3, 3))
    source_derivatives = np.empty((len(SOURCE_DEPTHS), len(ts), 3, 3, 3))

    for index, source_depth in enumerate(SOURCE_DEPTHS):
        G, dG_source, _ = pygrt.utils.lamb2(
            nu=NU, tbar=ts, R=EPICENTRAL_DISTANCE, depsrc=source_depth, azimuth=AZIMUTH
        )
        greens[index] = G
        source_derivatives[index] = dG_source

    return ts, greens, source_derivatives


def plot_components(
    ts: np.ndarray,
    data: np.ndarray,
    scale: float,
    offsets: np.ndarray,
    ylim: tuple[float, float],
    filename: str,
    derivative_coordinate: int | None = None,
) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(8.0, 7.0), layout="constrained")
    for axis, (i, j, component) in zip(axes.flat, COMPONENTS):
        for waveform, offset in zip(data, offsets):
            axis.plot(ts, scale * waveform[:, i, j] + offset, color="black", linewidth=0.8)
        axis.set_xlim(0.0, TMAX)
        axis.set_ylim(*ylim)
        axis.set_xticks(np.arange(0.0, TMAX + 0.01, 0.5))
        axis.tick_params(direction="in", top=True, right=True, length=3.0, width=0.7)
        axis.set_yticklabels([])
        axis.text(
            0.08,
            0.95,
            component_label(component, derivative_coordinate),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=11,
        )

    for source_depth, offset in zip(SOURCE_DEPTHS, offsets):
        axes[0, 0].text(
            0.03,
            offset,
            f"{source_depth:g}",
            transform=axes[0, 0].get_yaxis_transform(),
            ha="left",
            va="bottom",
            fontsize=9,
        )
    for axis in axes[-1, :]:
        axis.set_xlabel(r"$\bar{t}$", fontsize=11)

    fig.savefig(filename)
    plt.close(fig)


ts, greens, source_derivatives = compute_lamb2_results()
plot_components(ts, greens, GREEN_SCALE, GREEN_OFFSETS, GREEN_YLIM, "lamb2.svg")

for coordinate in (1, 2, 3):
    plot_components(
        ts,
        source_derivatives[:, :, coordinate - 1],
        DERIVATIVE_SCALES[coordinate],
        DERIVATIVE_OFFSETS,
        DERIVATIVE_YLIM,
        f"lamb2_d{coordinate}.svg",
        derivative_coordinate=coordinate,
    )
