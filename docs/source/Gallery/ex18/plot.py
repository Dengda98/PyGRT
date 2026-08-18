#!/usr/bin/env python3

from pathlib import Path
from typing import Dict, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf_file


HERE = Path(__file__).resolve().parent
CHANNELS = ("Z", "N", "E")
ROWS = ("PyGRT", "Okada", "PyGRT - Okada")


def read_nc(path: Path) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray], int]:
    """读取网格坐标、位移变量和位移偏导变量"""

    with netcdf_file(path, mode="r", mmap=False) as dataset:
        data = {name: np.array(variable[:], copy=True) for name, variable in dataset.variables.items()}
        rot2zne = int(getattr(dataset, "rot2ZNE", 0))
    if "north" not in data or "east" not in data:
        raise RuntimeError(f"{path} 不是网格形式的静态结果")
    return data["north"], data["east"], data, rot2zne


def check_results(
    north_a: np.ndarray,
    east_a: np.ndarray,
    north_b: np.ndarray,
    east_b: np.ndarray,
    data_a: Dict[str, np.ndarray],
    data_b: Dict[str, np.ndarray],
    names: Tuple[str, ...],
) -> None:
    """检查两个结果是否可以逐点比较"""

    if not np.allclose(north_a, north_b) or not np.allclose(east_a, east_b):
        raise RuntimeError("两个 NetCDF 文件的网格坐标不一致")
    for name in names:
        if name not in data_a or name not in data_b:
            raise RuntimeError(f"缺少变量 {name}")
        if data_a[name].shape != data_b[name].shape:
            raise RuntimeError(f"变量 {name} 的形状不一致")


def plot_compare(
    north: np.ndarray,
    east: np.ndarray,
    pygrt: Dict[str, np.ndarray],
    okada: Dict[str, np.ndarray],
    names: Tuple[str, ...],
    output: Path,
    title: str,
    colorbar_label: str,
) -> None:
    """绘制三行三列的结果对比图"""

    fig, axes = plt.subplots(3, 3, figsize=(9, 9), squeeze=False, layout="constrained")
    extent = (east.min(), east.max(), north.min(), north.max())
    images = []
    for column, name in enumerate(names):
        diff = pygrt[name] - okada[name]
        fields = (pygrt[name], okada[name], diff)
        limit = max(np.nanmax(np.abs(field)) for field in fields)
        if limit == 0.0:
            limit = 1.0
        for row, field in enumerate(fields):
            image = axes[row, column].imshow(
                field,
                extent=extent,
                origin="lower",
                aspect="equal",
                cmap="seismic",
                vmin=-limit,
                vmax=limit,
            )
            if row == 0:
                axes[row, column].set_title(name)
            if row == 2:
                axes[row, column].set_xlabel("East (km)")
            if column == 0:
                axes[row, column].set_ylabel(f"{ROWS[row]}\nNorth (km)")
        images.append(image)
        print(f"{name}: max difference = {np.max(np.abs(diff)):.6e}")
    for column in range(3):
        cax = axes[2, column].inset_axes([0.0, -0.25, 1.0, 0.05])
        colorbar = fig.colorbar(images[column], cax=cax, orientation="horizontal")
        colorbar.set_label(colorbar_label, fontsize=8, labelpad=2)
        colorbar.ax.tick_params(labelsize=7, pad=1)
    fig.suptitle(title)
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {output}")


derivative_names = tuple(f"{direction.lower()}{component}" for direction in CHANNELS for component in CHANNELS)

def compare_case(
    pygrt_path: Path,
    okada_path: Path,
    output_prefix: str,
    title_prefix: str,
) -> None:
    """检查并绘制一组 PyGRT 与 Okada 结果"""

    north_okada, east_okada, okada, okada_rot2zne = read_nc(okada_path)
    north_pygrt, east_pygrt, pygrt, pygrt_rot2zne = read_nc(pygrt_path)
    if okada_rot2zne != pygrt_rot2zne:
        raise RuntimeError(f"{title_prefix} 的输出坐标系不一致")
    if okada_rot2zne != 1:
        raise RuntimeError("本示例要求输出 ZNE 分量")

    check_results(
        north_okada,
        east_okada,
        north_pygrt,
        east_pygrt,
        okada,
        pygrt,
        (*CHANNELS, *derivative_names),
    )

    output_name = f"compare_{output_prefix}_displacement.svg" if output_prefix else "compare_displacement.svg"
    plot_compare(
        north_okada,
        east_okada,
        pygrt,
        okada,
        CHANNELS,
        HERE / output_name,
        f"{title_prefix} displacement comparison (ZNE)",
        "Displacement (cm)",
    )

    for direction in CHANNELS:
        names = tuple(f"{direction.lower()}{component}" for component in CHANNELS)
        output_name = f"compare_{output_prefix}_upar_{direction.lower()}.svg" if output_prefix else f"compare_upar_{direction.lower()}.svg"
        plot_compare(
            north_okada,
            east_okada,
            pygrt,
            okada,
            names,
            HERE / output_name,
            f"{title_prefix} displacement derivative comparison: d/d{direction.lower()}",
            "Dimensionless",
        )


compare_case(HERE / "static_syn.nc", HERE / "okada.nc", "", "Point-source")
compare_case(HERE / "finite_static_syn.nc", HERE / "finite_okada.nc", "finite", "Finite-fault")
