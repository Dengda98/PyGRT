#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from scipy.io import netcdf_file


def read_fault_surface(path: Path, fault_index: int) -> Tuple[np.ndarray, np.ndarray]:
    """读取一个有限接收断层的坐标和 ZNE 位移"""

    with netcdf_file(path, mmap=False) as dataset:
        norths = np.array(dataset.variables["north"].data, copy=True)
        easts = np.array(dataset.variables["east"].data, copy=True)
        depths = np.array(dataset.variables["depth"].data, copy=True)
        vertical_displacement = np.array(dataset.variables["Z"].data, copy=True)
        north_displacement = np.array(dataset.variables["N"].data, copy=True)
        east_displacement = np.array(dataset.variables["E"].data, copy=True)
        offsets = np.array(dataset.variables["offset"].data, copy=True)
        stksizes = np.array(dataset.variables["stksize"].data, copy=True)
        dipsizes = np.array(dataset.variables["dipsize"].data, copy=True)

    start = 0 if fault_index == 0 else int(offsets[fault_index - 1])
    end = int(offsets[fault_index])
    nstrike = int(stksizes[fault_index])
    ndip = int(dipsizes[fault_index])

    # 坐标使用 East、North、Depth 顺序，位移使用 Z、N、E 顺序
    coordinates = np.column_stack((easts[start:end], norths[start:end], depths[start:end]))
    displacements = np.column_stack((
        vertical_displacement[start:end],
        north_displacement[start:end],
        east_displacement[start:end],
    ))
    return (
        coordinates.reshape(ndip, nstrike, 3),
        displacements.reshape(ndip, nstrike, 3),
    )


def centers_to_edges(centers: np.ndarray) -> np.ndarray:
    """根据网格中心坐标生成 pcolormesh 所需的边界坐标"""

    centers = np.asarray(centers)
    if centers.size == 1:
        return np.array((0.0, 1.0))

    edges = np.empty(centers.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    # 将第一个网格边界平移到零点
    edges -= edges[0]
    return edges


def surface_axes(coordinates: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """根据相邻子断层中心的三维距离生成走向和倾向坐标"""

    strike_steps = np.linalg.norm(np.diff(coordinates[0, :, :], axis=0), axis=1)
    dip_steps = np.linalg.norm(np.diff(coordinates[:, 0, :], axis=0), axis=1)
    strike_centers = np.concatenate(([0.0], np.cumsum(strike_steps)))
    dip_centers = np.concatenate(([0.0], np.cumsum(dip_steps)))
    return centers_to_edges(strike_centers), centers_to_edges(dip_centers)


def plot_receiver_fault(
    input_path: Path,
    output_path: Path,
    fault_index: int = 0,
) -> None:
    """绘制一个有限接收断层的 ZNE 位移分量"""

    coordinates, displacements = read_fault_surface(Path(input_path), fault_index)
    strike_edges, dip_edges = surface_axes(coordinates)
    maximum = float(np.max(np.abs(displacements)))
    maximum = max(maximum, 1.0e-12)
    normal = Normalize(vmin=-maximum, vmax=maximum)
    colormap = plt.get_cmap("seismic")

    figure, axes = plt.subplots(1, 3, figsize=(13.5, 4.8))
    figure.subplots_adjust(left=0.06, right=0.98, bottom=0.22, top=0.88, wspace=0.25)
    component_names = ("Z", "N", "E")
    for icomponent, (axis, name) in enumerate(zip(axes, component_names)):
        # 三个分量使用同一个对称色标范围，便于比较正负位移
        mesh = axis.pcolormesh(
            strike_edges,
            dip_edges,
            displacements[:, :, icomponent],
            cmap=colormap,
            norm=normal,
            shading="flat",
        )
        axis.set_title(name)
        axis.set_xlabel("Along-strike distance (km)")
        if icomponent == 0:
            axis.set_ylabel("Along-dip depth (km)")
        axis.invert_yaxis()
        axis.set_aspect("equal", adjustable="box")

    # 将三个子图共用的水平色标放在图形中间下方
    colorbar_axis = figure.add_axes((0.31, 0.17, 0.38, 0.04))
    colorbar = figure.colorbar(mesh, cax=colorbar_axis, orientation="horizontal")
    colorbar.set_label("Displacement (cm)")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot finite receiver-fault ZNE displacements")
    parser.add_argument("input", nargs="?", type=Path, default=Path("stsyn_rf.nc"))
    parser.add_argument("output", nargs="?", type=Path, default=Path("rcv_fault_zne.svg"))
    parser.add_argument(
        "--fault-index",
        type=int,
        default=0,
        help="zero-based finite receiver-fault index to plot",
    )
    args = parser.parse_args()
    plot_receiver_fault(args.input, args.output, args.fault_index)


if __name__ == "__main__":
    main()
