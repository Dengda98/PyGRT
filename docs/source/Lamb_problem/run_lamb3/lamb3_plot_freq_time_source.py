"""比较 lamb3 源点导数和 greenfn 的 EX、DD、DS、SS 分量"""

import matplotlib

matplotlib.use("Agg")

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pygrt
from obspy import Stream, Trace, read

plt.rcParams.update({
    "font.sans-serif": "Times New Roman",
    "mathtext.fontset": "cm"
})

VP = 8.0  # km/s
VS = 4.62  # km/s
RHO = 3.3  # g/cm^3
NU = 0.5 * (1.0 - 2.0 * (VS / VP) ** 2) / (1.0 - (VS / VP) ** 2)
MU = VS**2 * RHO

SOURCE_DEPTH = 5.0  # km
RECEIVER_DEPTH = 0.1  # km
EPICENTRAL_DISTANCE = 10.0  # km
STRAIGHT_DISTANCE = np.hypot(EPICENTRAL_DISTANCE, SOURCE_DEPTH - RECEIVER_DEPTH)
NT = 1300
DT = 0.005  # s
PLOT_TMAX = 2.0
Y_LIMIT_MEAN_FACTOR = 4.0
Y_LIMIT_PERCENTILE = 98.0
Y_LIMIT_MARGIN = 1.1

MODEL_PATH = Path("_halfspace_mod")
GREEN_ROOT = Path("GRN")
OUTPUT_PATH = Path("lamb3_compare_freq_time_source.svg")

GREEN_CHANNELS = {"Z": "Z", "N": "R", "E": "T"}
SOURCE_COMPONENTS = {
    "EX": ("Z", "N"),
    "DD": ("Z", "N"),
    "DS": ("Z", "N", "E"),
    "SS": ("Z", "N", "E"),
}
SOURCE_NAMES = tuple(SOURCE_COMPONENTS)
OUTPUT_COMPONENTS = ("Z", "N", "E")


def remove_calculation_results() -> None:
    """删除脚本生成的中间结果"""
    if GREEN_ROOT.is_dir():
        shutil.rmtree(GREEN_ROOT)
    if MODEL_PATH.is_file():
        MODEL_PATH.unlink()


def convert_source_derivatives(dG_source: np.ndarray) -> np.ndarray:
    """将 lamb3 的 Gij,k' 转为 greenfn 的 EX、DD、DS、SS 分量

    lamb3 的坐标顺序为 x1=R、x2=T、x3=Z_down，数组索引为
    [时间, 源点导数方向, 接收点分量, 源点分量]
    """
    if dG_source.ndim != 4 or dG_source.shape[1:] != (3, 3, 3):
        raise ValueError("dG_source must have shape (nt, 3, 3, 3).")

    # 先在 lamb3 的 x1、x2、x3 坐标中构造四种震源分量
    nt = dG_source.shape[0]
    converted = np.zeros((nt, len(SOURCE_NAMES), 3))

    # EX = G_i1,1' + G_i2,2' + G_i3,3'
    converted[:, 0] = (
        dG_source[:, 0, :, 0] + dG_source[:, 1, :, 1] + dG_source[:, 2, :, 2]
    )

    # DD = 2 G_i3,3' - G_i1,1' - G_i2,2'
    converted[:, 1] = (
        2.0 * dG_source[:, 2, :, 2]
        - dG_source[:, 0, :, 0]
        - dG_source[:, 1, :, 1]
    )

    # DS 的 P-SV 部分使用 x1-z 剪切，SH 部分使用 x2-z 剪切
    ds_radial = -(
        dG_source[:, 0, :, 2] + dG_source[:, 2, :, 0]
    )
    ds_transverse = -(
        dG_source[:, 1, :, 2] + dG_source[:, 2, :, 1]
    )
    converted[:, 2, 0] = ds_radial[:, 0]
    converted[:, 2, 1] = ds_transverse[:, 1]
    converted[:, 2, 2] = ds_radial[:, 2]

    # SS 的 P-SV 部分为水平法向差，SH 部分为水平剪切
    ss_radial = dG_source[:, 0, :, 0] - dG_source[:, 1, :, 1]
    ss_transverse = dG_source[:, 0, :, 1] + dG_source[:, 1, :, 0]
    converted[:, 3, 0] = ss_radial[:, 0]
    converted[:, 3, 1] = ss_transverse[:, 1]
    converted[:, 3, 2] = ss_radial[:, 2]

    # 转为 greenfn 的 ZRT 输出顺序对应的 ZNE 顺序，并将 Z_down 改为 Z_up
    receiver_indices = (2, 0, 1)
    receiver_signs = np.array((-1.0, 1.0, 1.0))
    return converted[:, :, receiver_indices] * receiver_signs[None, None, :]


def calculate_lamb3() -> tuple[np.ndarray, np.ndarray]:
    """计算 lamb3 源点导数，并返回 tbar 和 ZNE 分量"""
    tbar = np.arange(NT, dtype=float) * DT * VS / STRAIGHT_DISTANCE
    _, dG_source, _, _ = pygrt.utils.lamb3(
        nu=NU,
        tbar=tbar,
        R=EPICENTRAL_DISTANCE,
        depsrc=SOURCE_DEPTH,
        deprcv=RECEIVER_DEPTH,
        azimuth=0.0,
    )
    return tbar, convert_source_derivatives(dG_source)


def calculate_greenfn() -> tuple[Stream, float]:
    """计算 greenfn 的四种源点导数基分量"""
    model = np.array([[0.0, VP, VS, RHO]])
    np.savetxt(MODEL_PATH, model)

    pymod = pygrt.PyModel1D(grn=GREEN_ROOT, modelpath=MODEL_PATH)
    pymod.greenfn(
        depsrc=SOURCE_DEPTH,
        deprcv=RECEIVER_DEPTH,
        dists=[EPICENTRAL_DISTANCE],
        nt=NT,
        dt=DT,
        Length=50.0,
        gf_source=["EX", "DC"],
        print_log=False,
    )
    st = read(str(GREEN_ROOT / "*" / "*.sac"))

    # greenfn 输出的是脉冲型格林函数，积分一次后与 lamb3 的阶跃型结果对应
    pygrt.utils.stream_integral(st)

    scale = np.pi**2 * MU * STRAIGHT_DISTANCE**2
    return st, scale


def select_trace(st: Stream, source: str, component: str) -> Trace:
    """读取一个 greenfn 源类型和接收分量"""
    channel = source + GREEN_CHANNELS[component]
    traces = st.select(channel=channel)
    if len(traces) != 1:
        raise ValueError(f"Expected one {channel} trace, got {len(traces)}.")
    return traces[0]


def determine_y_limits(*series: np.ndarray) -> tuple[float, float]:
    """根据序列的典型振幅确定对称的 y 轴范围"""
    values = np.concatenate([np.asarray(item).ravel() for item in series])
    values = values[np.isfinite(values)]
    if values.size == 0:
        return -1.0, 1.0

    absolute_values = np.abs(values)
    mean_absolute_value = np.mean(absolute_values)
    percentile_value = np.percentile(absolute_values, Y_LIMIT_PERCENTILE)
    typical_amplitude = max(Y_LIMIT_MEAN_FACTOR * mean_absolute_value, percentile_value)
    typical_amplitude = max(typical_amplitude, np.finfo(float).eps)
    limit = Y_LIMIT_MARGIN * typical_amplitude
    return -limit, limit


def plot_comparison(tbar: np.ndarray, lamb3_data: np.ndarray, st: Stream, scale: float) -> None:
    """绘制 5×4 的源点导数对比图"""
    fig, axes = plt.subplots(5, 4, figsize=(13.0, 10.0), sharex=True)
    fig.subplots_adjust(
        left=0.06,
        right=0.96,
        bottom=0.07,
        top=0.94,
        hspace=0.15,
        wspace=0.08,
    )

    # 加大左右两组子图之间的间隔
    middle_gap = 0.02
    for ax in axes[:, 2:].flat:
        position = ax.get_position()
        ax.set_position(
            [position.x0 + middle_gap, position.y0, position.width, position.height]
        )

    plot_rows = (
        (("EX", "Z"), ("DS", "N")),
        (("EX", "N"), ("DS", "E")),
        (("DD", "Z"), ("SS", "Z")),
        (("DD", "N"), ("SS", "N")),
        (("DS", "Z"), ("SS", "E")),
    )
    plot_mask = (tbar >= 0.0) & (tbar <= PLOT_TMAX)

    source_indices = {source: index for index, source in enumerate(SOURCE_NAMES)}
    component_indices = {
        component: index for index, component in enumerate(OUTPUT_COMPONENTS)
    }

    for row, (left_item, right_item) in enumerate(plot_rows):
        for group, (source, component) in enumerate((left_item, right_item)):
            source_index = source_indices[source]
            component_index = component_indices[component]
            lamb3_axis = axes[row, 2 * group]
            comparison_axis = axes[row, 2 * group + 1]
            trace = select_trace(st, source, component)
            greenfn_data = trace.data * scale
            lamb3_series = lamb3_data[:, source_index, component_index]
            y_limits = determine_y_limits(lamb3_series[plot_mask], greenfn_data[plot_mask])

            lamb3_axis.plot(tbar, lamb3_series, color="0.5", linewidth=1.2)
            comparison_axis.plot(tbar, lamb3_series, color="0.5", linewidth=1.2)
            comparison_axis.plot(tbar, greenfn_data, color="blue", linewidth=0.7)
            for ax in (lamb3_axis, comparison_axis):
                ax.set_xlim(0.0, PLOT_TMAX)
                ax.set_ylim(*y_limits)
                ax.grid(linewidth=0.4)
                ax.tick_params(direction="in", top=True, right=True, length=3.0, width=0.7)
                ax.text(
                    0.05,
                    0.92,
                    f"{source}{GREEN_CHANNELS[component]}",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=11,
                )

            comparison_axis.tick_params(labelleft=False)

    for column in range(4):
        axes[0, column].set_title(
            "From Time-Domain" if column % 2 == 0 else "From Frequency-Domain", fontsize=12
        )
        axes[-1, column].set_xlabel(r"$\bar{t}$", fontsize=11)

    fig.savefig(OUTPUT_PATH, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """执行计算、绘图并清理中间文件"""
    remove_calculation_results()
    try:
        tbar, lamb3_data = calculate_lamb3()
        st, scale = calculate_greenfn()
        plot_comparison(tbar, lamb3_data, st, scale)
    finally:
        remove_calculation_results()


if __name__ == "__main__":
    main()
