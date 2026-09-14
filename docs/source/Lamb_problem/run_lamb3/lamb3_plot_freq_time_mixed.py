"""比较 lamb3 二阶混合导数和 greenfn 的 z/r 前缀分量"""

import matplotlib

matplotlib.use("Agg")

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
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
FREQUENCY_BASELINE_TMAX = 0.52
INTEGRATION_COUNT = 1   # 关键参数，用于覆盖频域解的 Gibbs 效应
Y_LIMIT_MEAN_FACTOR = 4.0
Y_LIMIT_PERCENTILE = 98.0
Y_LIMIT_MARGIN = 1.1

P_WAVE_ARRIVAL = np.sqrt(0.5 * (1.0 - 2.0 * NU) / (1.0 - NU))
S_WAVE_ARRIVAL = 1.0
# 只按直达 P、S 两个主到时分段，拟合窗口避开其他震相波前
FREQUENCY_ALIGNMENT_SEGMENTS = (
    (P_WAVE_ARRIVAL, S_WAVE_ARRIVAL, P_WAVE_ARRIVAL + 0.07, 0.85),
    (S_WAVE_ARRIVAL, np.inf, S_WAVE_ARRIVAL + 0.10, 1.45),
)

MODEL_PATH = Path("_halfspace_mod")
GREEN_ROOT = Path("GRN")
OUTPUT_PATHS = {
    "z": Path("lamb3_compare_freq_time_mixed_z.svg"),
    "r": Path("lamb3_compare_freq_time_mixed_r.svg"),
}

GREEN_CHANNELS = {"Z": "Z", "N": "R", "E": "T"}
SOURCE_COMPONENTS = {
    "EX": ("Z", "N"),
    "DD": ("Z", "N"),
    "DS": ("Z", "N", "E"),
    "SS": ("Z", "N", "E"),
}
SOURCE_NAMES = tuple(SOURCE_COMPONENTS)
OUTPUT_COMPONENTS = ("Z", "N", "E")
AXIS_MAP = (2, 0, 1)
AXIS_SIGN = np.array((-1.0, 1.0, 1.0))


def remove_calculation_results() -> None:
    """删除脚本生成的中间结果"""
    if GREEN_ROOT.is_dir():
        shutil.rmtree(GREEN_ROOT)
    if MODEL_PATH.is_file():
        MODEL_PATH.unlink()


def convert_mixed_derivatives(dG_mixed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """将 lamb3 的混合导数转为 greenfn 的 z/r EX、DD、DS、SS 分量

    lamb3 的坐标顺序为 x1=R、x2=T、x3=Z_down，方位角为 0 时对应
    全局坐标 N、E、Z_down，数组索引为
    [时间, 接收点导数方向, 源点导数方向, 接收点分量, 源点分量]
    """
    if dG_mixed.ndim != 5 or dG_mixed.shape[1:] != (3, 3, 3, 3):
        raise ValueError("dG_mixed must have shape (nt, 3, 3, 3, 3).")

    # 将接收点导数、源点导数、接收分量和源分量全部转换为 ZNE
    mapped = dG_mixed
    for axis in (1, 2, 3, 4):
        mapped = np.take(mapped, AXIS_MAP, axis=axis)
    signs = (
        AXIS_SIGN[None, :, None, None, None]
        * AXIS_SIGN[None, None, :, None, None]
        * AXIS_SIGN[None, None, None, :, None]
        * AXIS_SIGN[None, None, None, None, :]
    )
    mapped = mapped * signs

    # 在 ZNE 坐标中构造四种 greenfn 震源分量
    nt = dG_mixed.shape[0]
    converted = np.zeros((nt, len(SOURCE_NAMES), 3, 3))

    # EX = G_i1,1' + G_i2,2' + G_i3,3'
    converted[:, 0] = (
        mapped[:, :, 0, :, 0]
        + mapped[:, :, 1, :, 1]
        + mapped[:, :, 2, :, 2]
    )

    # DD = 2 G_i3,3' - G_i1,1' - G_i2,2'
    converted[:, 1] = (
        2.0 * mapped[:, :, 0, :, 0]
        - mapped[:, :, 1, :, 1]
        - mapped[:, :, 2, :, 2]
    )

    # DS 的 P-SV 部分使用 N-Z 剪切，SH 部分使用 E-Z 剪切
    ds_n = mapped[:, :, 0, :, 1] + mapped[:, :, 1, :, 0]
    ds_e = mapped[:, :, 0, :, 2] + mapped[:, :, 2, :, 0]
    converted[:, 2, :, 0] = ds_n[:, :, 0]
    converted[:, 2, :, 1] = ds_n[:, :, 1]
    converted[:, 2, :, 2] = ds_e[:, :, 2]

    # SS 的 Z/R 通道沿用 greenfn 的水平法向差约定，T 通道为水平剪切
    ss_n = mapped[:, :, 1, :, 1] - mapped[:, :, 2, :, 2]
    ss_e = mapped[:, :, 2, :, 1] + mapped[:, :, 1, :, 2]
    converted[:, 3, :, 0] = ss_n[:, :, 0]
    converted[:, 3, :, 1] = ss_n[:, :, 1]
    converted[:, 3, :, 2] = ss_e[:, :, 2]

    # ZNE 中 0、1 两个接收点导数方向分别对应 greenfn 的 z、r 前缀
    z_data = converted[:, :, 0]
    r_data = converted[:, :, 1]
    return z_data, r_data


def integrate_time_series(data: np.ndarray, step: float, count: int) -> np.ndarray:
    """按梯形法对时间序列积分指定次数"""
    result = np.array(data, dtype=float, copy=True)
    for _ in range(count):
        last = result[0].copy()
        result[0] = 0.0
        for index in range(1, len(result)):
            current = result[index].copy()
            result[index] = 0.5 * (current + last) * step + result[index - 1]
            last = current
    return result


def remove_frequency_baseline(st: Stream, tbar: np.ndarray) -> None:
    """去掉二次频域积分在首个 P 波前产生的一次基线漂移"""
    mask = tbar < FREQUENCY_BASELINE_TMAX
    for trace in st:
        coefficients = np.polyfit(tbar[mask], trace.data[mask], INTEGRATION_COUNT)
        trace.data[:] -= np.polyval(coefficients, tbar)


def align_frequency_series(tbar: np.ndarray, reference: np.ndarray, actual: np.ndarray) -> np.ndarray:
    """按直达 P、S 两个主到时对齐频域积分后的低频漂移"""
    aligned = np.array(actual, dtype=float, copy=True)
    for left, right, fit_left, fit_right in FREQUENCY_ALIGNMENT_SEGMENTS:
        apply_mask = (tbar >= left) & (tbar < right)
        fit_mask = apply_mask & (tbar >= fit_left) & (tbar <= fit_right)
        if np.count_nonzero(fit_mask) < 3:
            continue

        design = np.stack((np.ones(np.count_nonzero(fit_mask)), tbar[fit_mask]), axis=1)
        coefficients = np.linalg.lstsq(
            design,
            reference[fit_mask] - actual[fit_mask],
            rcond=None,
        )[0]
        apply_design = np.stack((np.ones(np.count_nonzero(apply_mask)), tbar[apply_mask]), axis=1)
        aligned[apply_mask] += apply_design @ coefficients
    return aligned


def calculate_lamb3() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """计算 lamb3 混合导数，并返回 tbar 和 ZNE 分量"""
    tbar = np.arange(NT, dtype=float) * DT * VS / STRAIGHT_DISTANCE
    _, _, _, dG_mixed = pygrt.utils.lamb3(
        nu=NU,
        tbar=tbar,
        R=EPICENTRAL_DISTANCE,
        depsrc=SOURCE_DEPTH,
        deprcv=RECEIVER_DEPTH,
        azimuth=0.0,
    )
    z_data, r_data = convert_mixed_derivatives(dG_mixed)
    step = tbar[1] - tbar[0]
    z_data = integrate_time_series(z_data, step, INTEGRATION_COUNT)
    r_data = integrate_time_series(r_data, step, INTEGRATION_COUNT)
    return tbar, z_data, r_data


def calculate_greenfn() -> tuple[Stream, float]:
    """计算 greenfn 的 z/r 二阶混合导数基分量"""
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
        calc_upar=True,
        gf_source=["EX", "DC"],
        print_log=False,
    )
    st = read(str(GREEN_ROOT / "*" / "*.sac"))

    # 第一次积分得到阶跃响应，额外一次积分得到与 lamb3 一次积分结果对应的量
    for _ in range(INTEGRATION_COUNT + 1):
        pygrt.utils.stream_integral(st)

    frequency_tbar = np.arange(NT, dtype=float) * DT * VS / STRAIGHT_DISTANCE
    remove_frequency_baseline(st, frequency_tbar)
    scale = np.pi**2 * MU * STRAIGHT_DISTANCE**3 * (VS / STRAIGHT_DISTANCE)**INTEGRATION_COUNT
    return st, scale


def select_trace(st: Stream, prefix: str, source: str, component: str) -> Trace:
    """读取一个 greenfn 导数前缀、源类型和接收分量"""
    channel = f"{prefix}{source}{GREEN_CHANNELS[component]}"
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


def plot_comparison(
    tbar: np.ndarray,
    lamb3_data: np.ndarray,
    st: Stream,
    scale: float,
    prefix: str,
) -> None:
    """绘制 5×4 的混合导数对比图"""
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
            trace = select_trace(st, prefix, source, component)
            lamb3_series = lamb3_data[:, source_index, component_index]
            greenfn_data = align_frequency_series(tbar, lamb3_series, trace.data * scale)
            y_limits = determine_y_limits(lamb3_series[plot_mask], greenfn_data[plot_mask])

            for ax in (lamb3_axis, comparison_axis):
                ax.plot(tbar, lamb3_series, color="0.5", linewidth=1.2)
                ax.set_xlim(0.0, PLOT_TMAX)
                ax.set_ylim(*y_limits)
                ax.grid(linewidth=0.4)
                ax.tick_params(direction="in", top=True, right=True, length=3.0, width=0.7)
                ax.text(
                    0.05,
                    0.92,
                    f"{prefix}{source}{GREEN_CHANNELS[component]}",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=11,
                )

            comparison_axis.plot(tbar, greenfn_data, color="blue", linewidth=0.7)
            comparison_axis.tick_params(labelleft=False)

    for column in range(4):
        axes[0, column].set_title(
            "From Time-Domain" if column % 2 == 0 else "From Frequency-Domain (P/S aligned)", fontsize=12
        )
        axes[-1, column].set_xlabel(r"$\bar{t}$", fontsize=11)

    fig.savefig(OUTPUT_PATHS[prefix], bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """执行计算、绘图并清理中间文件"""
    remove_calculation_results()
    try:
        tbar, z_data, r_data = calculate_lamb3()
        st, scale = calculate_greenfn()
        plot_comparison(tbar, z_data, st, scale, "z")
        plot_comparison(tbar, r_data, st, scale, "r")
    finally:
        remove_calculation_results()


if __name__ == "__main__":
    main()
