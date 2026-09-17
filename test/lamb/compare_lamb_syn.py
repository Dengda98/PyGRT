"""比较解析 Lamb 模块与 greenfn+syn 模块的物理位移及空间偏导结果"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
from obspy import read


VP = 8.0
VS = 4.62
RHO = 3.3
NU = 0.5 * (1.0 - 2.0 * (VS / VP) ** 2) / (1.0 - (VS / VP) ** 2)
R = 10.0
AZIMUTH = 30.0
NT = 1300
DT = 0.005
# 是否绘制对比结果
DOPLOT = False
PLOT_DIRECTORY = Path(__file__).resolve().parent / "compare_lamb_syn_plots"
# 对偏导记录额外积分一次，以减弱脉冲尖峰对对比和绘图的影响
DERIVATIVE_EXTRA_INTEGRATIONS = 1
# 对振幅和误差的 L2 范数进行高分位截断，避免孤立尖峰主导评分
ROBUST_PERCENTILE = 95.0
MAX_ROBUST_RELATIVE_L2_ERROR = 0.1
P_WAVE_ARRIVAL = np.sqrt(0.5 * (1.0 - 2.0 * NU) / (1.0 - NU))
S_WAVE_ARRIVAL = 1.0
# 根据 syn 的积分阶数和额外积分，低频漂移在偏导记录中表现为相应阶数基线
FREQUENCY_BASELINE_TMAX = P_WAVE_ARRIVAL - 0.05
ALIGNMENT_SEGMENTS = (
    (P_WAVE_ARRIVAL, S_WAVE_ARRIVAL, P_WAVE_ARRIVAL + 0.07, 0.85),
    (S_WAVE_ARRIVAL, np.inf, S_WAVE_ARRIVAL + 0.10, 1.45),
)


@dataclass(frozen=True)
class SourceCase:
    name: str
    options: tuple[str, ...]
    scale: str
    is_moment: bool

    @property
    def integration_order(self) -> int:
        # 单力源使用一次积分，偶极和力偶源使用两次积分
        return 2 if self.is_moment else 1


SF_SOURCE_CASE = SourceCase("SF", ("-F0.8/-0.4/1.2",), "1e15", False)
SOURCE_CASES = (
    SourceCase("EX", (), "1e20", True),
    SF_SOURCE_CASE,
    SourceCase("DC", ("-M100/30/70",), "1e20", True),
    SourceCase("TS", ("-M100/30",), "1e20", True),
    SourceCase("MT", ("-T1/-2/3/0.5/1.2/-0.7",), "1e20", True),
)


@dataclass(frozen=True)
class GeometryCase:
    name: str
    source_depth: float
    receiver_depth: float
    lamb_depth_options: tuple[str, ...]
    greenfn_depth: str
    source_cases: tuple[SourceCase, ...] = SOURCE_CASES
    validate_derivatives: bool = True


GEOMETRY_CASES = (
    GeometryCase("lamb2", 5.0, 0.0, ("-Ds5", "-Dr0"), "5/0"),
    GeometryCase("lamb3", 5.0, 1.0, ("-Ds5", "-Dr1"), "5/1"),
    GeometryCase(
        "surface", 0.0, 0.0, ("-Ds0", "-Dr0"), "0/0",
        source_cases=(SF_SOURCE_CASE,), validate_derivatives=False,
    ),
)


def run(command: list[str]) -> None:
    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)


def find_greenfn_directory(root: Path) -> Path:
    candidates = [path for path in root.iterdir() if path.is_dir() and (path / "EXZ.sac").is_file()]
    if len(candidates) != 1:
        raise AssertionError(f"Expected one Green-function directory under {root}, got {candidates}")
    return candidates[0]


def integrate_series(data: np.ndarray, step: float, count: int) -> np.ndarray:
    """使用梯形法对时间序列积分指定次数"""
    if count < 0:
        raise ValueError("The integration count must be nonnegative")
    result = np.array(data, dtype=float, copy=True)
    for _ in range(count):
        integrated = np.empty_like(result)
        integrated[0] = 0.0
        integrated[1:] = np.cumsum(0.5 * (result[1:] + result[:-1]) * step)
        result = integrated
    return result


def read_trace(path: Path, integration_order: int = 0) -> np.ndarray:
    """读取 SAC 记录并在读取后按需积分"""
    if not path.is_file():
        raise AssertionError(f"Missing SAC file: {path}")
    trace = read(str(path))[0]
    if trace.stats.npts != NT or not np.isclose(trace.stats.delta, DT):
        raise AssertionError(f"Unexpected SAC header in {path}")
    data = np.asarray(trace.data, dtype=float)
    if not np.isfinite(data).all():
        raise AssertionError(f"Non-finite data in {path}")
    return integrate_series(data, trace.stats.delta, integration_order)


def remove_frequency_baseline(data: np.ndarray, tbar: np.ndarray, polynomial_order: int) -> np.ndarray:
    """使用首个 P 波前的静音段去除频域结果的低频基线"""
    if polynomial_order < 0:
        raise ValueError("The baseline polynomial order must be nonnegative")
    mask = tbar < FREQUENCY_BASELINE_TMAX
    if np.count_nonzero(mask) <= polynomial_order:
        raise ValueError("Not enough pre-P samples for frequency baseline removal")
    coefficients = np.polyfit(tbar[mask], data[mask], polynomial_order)
    return data - np.polyval(coefficients, tbar)


def align_frequency_series(lamb_data: np.ndarray, syn_data: np.ndarray, tbar: np.ndarray) -> np.ndarray:
    aligned = np.array(syn_data, dtype=float, copy=True)
    for left, right, fit_left, fit_right in ALIGNMENT_SEGMENTS:
        apply_mask = (tbar >= left) & (tbar < right)
        fit_mask = apply_mask & (tbar >= fit_left) & (tbar <= fit_right)
        if np.count_nonzero(fit_mask) < 3:
            continue
        design = np.stack((np.ones(np.count_nonzero(fit_mask)), tbar[fit_mask]), axis=1)
        coefficients = np.linalg.lstsq(
            design,
            lamb_data[fit_mask] - syn_data[fit_mask],
            rcond=None,
        )[0]
        apply_design = np.stack((np.ones(np.count_nonzero(apply_mask)), tbar[apply_mask]), axis=1)
        aligned[apply_mask] += apply_design @ coefficients
    return aligned


def robust_l2_norm(data: np.ndarray) -> float:
    """计算对孤立高峰不敏感的截断 L2 范数"""
    absolute = np.abs(np.asarray(data, dtype=float))
    if absolute.size == 0:
        return 0.0
    cap = np.percentile(absolute, ROBUST_PERCENTILE)
    return float(np.sqrt(np.mean(np.minimum(absolute, cap) ** 2)))


def compare_series_data(
    lamb_data: np.ndarray, syn_data: np.ndarray, time: np.ndarray, distance: float
) -> tuple[np.ndarray, float]:
    tbar = time * VS / distance
    aligned_syn_data = align_frequency_series(lamb_data, syn_data, tbar)
    compare_mask = (tbar >= P_WAVE_ARRIVAL + 0.08) & (np.abs(tbar - S_WAVE_ARRIVAL) > 0.08)
    if np.count_nonzero(compare_mask) < 10:
        raise AssertionError("The comparison window is too short")
    lamb_scale = robust_l2_norm(lamb_data[compare_mask])
    syn_scale = robust_l2_norm(aligned_syn_data[compare_mask])
    scale = np.sqrt(0.5 * (lamb_scale**2 + syn_scale**2))
    if scale <= np.finfo(float).eps:
        return aligned_syn_data, 0.0
    error = robust_l2_norm((lamb_data - aligned_syn_data)[compare_mask])
    return aligned_syn_data, error / scale


def compare_series(
    lamb_data: np.ndarray, syn_data: np.ndarray, time: np.ndarray, distance: float
) -> float:
    """计算两道记录的稳健归一化 L2 误差"""
    _, score = compare_series_data(lamb_data, syn_data, time, distance)
    return score


@dataclass(frozen=True)
class ComparisonResult:
    component: str
    lamb_data: np.ndarray
    syn_data: np.ndarray
    score: float


def determine_plot_limits(*series: np.ndarray) -> tuple[float, float]:
    """根据典型振幅设置对称的绘图范围"""
    values = np.concatenate([np.asarray(item, dtype=float).ravel() for item in series])
    values = values[np.isfinite(values)]
    if values.size == 0:
        return -1.0, 1.0
    typical_amplitude = np.percentile(np.abs(values), 99.0)
    limit = 1.1 * max(typical_amplitude, np.finfo(float).eps)
    return -limit, limit


def plot_comparison_case(
    geometry: GeometryCase, source: SourceCase, zne: bool,
    time: np.ndarray, comparisons: list[ComparisonResult],
) -> Path:
    """为一个 case 生成每个分量一页的多页 PDF"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    coordinate_name = "ZNE" if zne else "ZRT"
    PLOT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    pdf_path = PLOT_DIRECTORY / f"{geometry.name}_{source.name}_{coordinate_name}.pdf"
    with PdfPages(pdf_path) as pdf:
        for comparison in comparisons:
            figure, axis = plt.subplots(figsize=(10.0, 4.8))
            axis.plot(
                time, comparison.lamb_data, color="tab:blue", linewidth=1.4,
                linestyle="-", label="Lamb",
            )
            axis.plot(
                time, comparison.syn_data, color="tab:orange", linewidth=0.9,
                linestyle="--", label="greenfn + syn",
            )
            axis.set_xlim(time[0], time[-1])
            axis.set_ylim(*determine_plot_limits(comparison.lamb_data, comparison.syn_data))
            axis.set_xlabel("Time (s)")
            axis.set_ylabel("Amplitude")
            axis.set_title(
                f"{geometry.name} {source.name} {coordinate_name} {comparison.component} "
                f"(robust relative L2={comparison.score:.3f})"
            )
            axis.grid(linewidth=0.5)
            axis.legend()
            figure.tight_layout()
            pdf.savefig(figure)
            plt.close(figure)
    return pdf_path


def lamb_command(geometry: GeometryCase, source: SourceCase, output: Path, zne: bool) -> list[str]:
    command = [
        "grt",
        "lamb",
        f"-H{VP:.16g}/{VS:.16g}/{RHO:.16g}",
        f"-N{NT}/{DT}",
        f"-R{R}",
        f"-A{AZIMUTH}",
        *geometry.lamb_depth_options,
        f"-S{source.scale}",
        *source.options,
        f"-I{source.integration_order}",
        "-e",
        f"-O{output}",
        "-s",
    ]
    if zne:
        command.append("-n")
    return command


def syn_command(green_directory: Path, source: SourceCase, output: Path, zne: bool) -> list[str]:
    command = [
        "grt",
        "syn",
        f"-G{green_directory}",
        f"-A{AZIMUTH}",
        f"-S{source.scale}",
        *source.options,
        f"-I{source.integration_order}",
        "-e",
        f"-O{output}",
        "-s",
    ]
    if zne:
        command.append("-N")
    return command


def compare_case(geometry: GeometryCase, source: SourceCase, green_directory: Path, root: Path, zne: bool) -> list[tuple[str, float]]:
    coordinate_name = "ZNE" if zne else "ZRT"
    lamb_directory = root / f"lamb_{geometry.name}_{source.name}_{coordinate_name}"
    syn_directory = root / f"syn_{geometry.name}_{source.name}_{coordinate_name}"
    run(lamb_command(geometry, source, lamb_directory, zne))
    run(syn_command(green_directory, source, syn_directory, zne))

    channels = "ZNE" if zne else "ZRT"
    directions = "zne" if zne else "zrt"
    time = np.arange(NT, dtype=float) * DT
    distance = np.hypot(R, geometry.source_depth - geometry.receiver_depth)
    tbar = time * VS / distance
    scores = []
    comparisons = []
    for channel in channels:
        lamb_data = read_trace(lamb_directory / f"{channel}.sac")
        syn_data = read_trace(syn_directory / f"{channel}.sac")
        aligned_syn_data, score = compare_series_data(lamb_data, syn_data, time, distance)
        scores.append((channel, score))
        comparisons.append(ComparisonResult(channel, lamb_data, aligned_syn_data, score))
        if score > MAX_ROBUST_RELATIVE_L2_ERROR:
            raise AssertionError(
                f"{geometry.name} {source.name} {coordinate_name} {channel} error={score:.3f} "
                f"> threshold={MAX_ROBUST_RELATIVE_L2_ERROR:.3f}"
            )

    if geometry.validate_derivatives:
        # 对非地表情景，验证三个导数方向与三个接收分量的全部组合
        for direction in directions:
            for channel in channels:
                component = f"{direction}{channel}"
                lamb_data = read_trace(
                    lamb_directory / f"{component}.sac",
                    integration_order=DERIVATIVE_EXTRA_INTEGRATIONS,
                )
                syn_data = read_trace(
                    syn_directory / f"{component}.sac",
                    integration_order=DERIVATIVE_EXTRA_INTEGRATIONS,
                )
                # syn 的积分使漂移阶数增加，额外积分后基线阶数为总积分阶数减一
                baseline_order = source.integration_order + DERIVATIVE_EXTRA_INTEGRATIONS - 1
                syn_data = remove_frequency_baseline(syn_data, tbar, baseline_order)
                aligned_syn_data, score = compare_series_data(lamb_data, syn_data, time, distance)
                scores.append((component, score))
                comparisons.append(ComparisonResult(component, lamb_data, aligned_syn_data, score))
                if score > MAX_ROBUST_RELATIVE_L2_ERROR:
                    raise AssertionError(
                        f"{geometry.name} {source.name} {coordinate_name} {component} error={score:.3f} "
                        f"> threshold={MAX_ROBUST_RELATIVE_L2_ERROR:.3f}"
                    )
    else:
        # 地表—地表时 -e 被 Lamb 模块忽略，因此只允许输出三条位移记录
        expected_files = {f"{channel}.sac" for channel in channels}
        actual_files = {path.name for path in lamb_directory.glob("*.sac")}
        if actual_files != expected_files:
            raise AssertionError(
                f"{geometry.name} {source.name} should only write displacement files; "
                f"got {sorted(actual_files)}"
            )
    if DOPLOT:
        plot_comparison_case(geometry, source, zne, time, comparisons)
    return scores


def main() -> None:
    with TemporaryDirectory(prefix="pygrt-lamb-compare-") as temporary:
        root = Path(temporary)
        model = root / "halfspace"
        model.write_text(f"0.0 {VP} {VS} {RHO}\n", encoding="ascii")

        for geometry in GEOMETRY_CASES:
            green_root = root / f"green_{geometry.name}"
            run([
                "grt",
                "greenfn",
                f"-M{model}",
                f"-D{geometry.greenfn_depth}",
                f"-N{NT}/{DT}",
                f"-R{R}",
                "-e",
                f"-O{green_root}",
                "-s",
            ])
            green_directory = find_greenfn_directory(green_root)
            for source in geometry.source_cases:
                for zne in (False, True):
                    scores = compare_case(geometry, source, green_directory, root, zne)
                    worst_score = max(score for _, score in scores)
                    checked = "displacement" if not geometry.validate_derivatives else "displacement+derivatives"
                    print(
                        f"{geometry.name} {source.name} {'ZNE' if zne else 'ZRT'} "
                        f"checked={checked} worst_robust_rel_l2={worst_score:.3f}"
                    )

    print(
        "compare_lamb_syn.py: all comparisons passed, "
        f"threshold={MAX_ROBUST_RELATIVE_L2_ERROR:.3f}"
    )


if __name__ == "__main__":
    main()
