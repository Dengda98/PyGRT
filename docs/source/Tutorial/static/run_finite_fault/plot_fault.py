#!/usr/bin/env python3

from pathlib import Path
from typing import List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


Fault = Tuple[float, float, float, float, float, float, float]


def read_faults(path: Path) -> List[Fault]:
    """读取 Coulomb 格式断层文件中的断层几何"""

    faults = []
    with path.open(encoding="utf-8") as file:
        lines = file.readlines()

    for line_number, line in enumerate(lines[2:], start=3):
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) < 11:
            raise ValueError(f"{path}:{line_number}: expected 11 columns")

        try:
            values = [float(field) for field in fields[:11]]
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid numeric value") from exc

        _, east_begin, north_begin, east_end, north_end, _, _, _, dip, top, bot = values
        if dip <= 0.0 or dip > 90.0:
            raise ValueError(f"{path}:{line_number}: dip must be in (0, 90]")
        if bot <= top:
            raise ValueError(f"{path}:{line_number}: bot must be greater than top")

        faults.append((east_begin, north_begin, east_end, north_end, dip, top, bot))

    if not faults:
        raise ValueError(f"{path}: no fault geometry found")
    return faults


def fault_geometry(fault: Fault) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """计算断层顶边、底边及其三维顶点"""

    east_begin, north_begin, east_end, north_end, dip, top, bot = fault
    delta_east = east_end - east_begin
    delta_north = north_end - north_begin
    length = np.hypot(delta_east, delta_north)
    if length <= 0.0:
        raise ValueError("fault along-strike length must be positive")

    dip_rad = np.deg2rad(dip)
    horizontal_width = (bot - top) * np.cos(dip_rad) / np.sin(dip_rad)
    down_dip_shift = np.array([
        horizontal_width * delta_north / length,
        -horizontal_width * delta_east / length,
    ])

    top_edge = np.array([
        [east_begin, north_begin],
        [east_end, north_end],
    ])
    bottom_edge = top_edge + down_dip_shift
    projected = np.array([
        top_edge[0], top_edge[1], bottom_edge[1], bottom_edge[0],
    ])
    depths = np.array([top, top, bot, bot])
    vertices_3d = np.column_stack((projected, depths))
    return top_edge, bottom_edge, vertices_3d


def set_horizontal_limits(ax, points: np.ndarray) -> None:
    """设置水平投影坐标范围"""

    minimum = points.min(axis=0)
    maximum = points.max(axis=0)
    span = np.maximum(maximum - minimum, 1.0)
    padding = 0.12 * max(span)
    ax.set_xlim(minimum[0] - padding, maximum[0] + padding)
    ax.set_ylim(minimum[1] - padding, maximum[1] + padding)


def plot_fault_geometry(input_path: Path, output_path: Path) -> None:
    """绘制有限断层的水平投影和三维展布"""

    input_path = Path(input_path)
    output_path = Path(output_path)

    faults = read_faults(input_path)
    geometries = [fault_geometry(fault) for fault in faults]
    projected = np.vstack([
        np.vstack((top_edge, bottom_edge))
        for top_edge, bottom_edge, _ in geometries
    ])

    top_color = "#0072B2"
    bottom_color = "#D55E00"
    face_color = "#56B4E9"
    fig = plt.figure(figsize=(12, 5.2), constrained_layout=True)
    ax_map = fig.add_subplot(1, 2, 1)
    ax_3d = fig.add_subplot(1, 2, 2, projection="3d")

    top_handles = []
    bottom_handles = []
    for index, (fault, geometry) in enumerate(zip(faults, geometries), start=1):
        top_edge, bottom_edge, vertices_3d = geometry
        _, _, _, _, _, top, bot = fault
        polygon = vertices_3d[:, :2]

        ax_map.fill(
            polygon[:, 0], polygon[:, 1], color=face_color, alpha=0.18,
            zorder=1,
        )
        for top_point, bottom_point in zip(top_edge, bottom_edge):
            ax_map.plot(
                [top_point[0], bottom_point[0]],
                [top_point[1], bottom_point[1]],
                color="0.45", linestyle="--", linewidth=0.9, zorder=2,
            )
        top_handle, = ax_map.plot(
            top_edge[:, 0], top_edge[:, 1], "o-", color=top_color,
            linewidth=2.0, markersize=5.5,
            label=f"Fault {index} top ({top:g} km)", zorder=3,
        )
        bottom_handle, = ax_map.plot(
            bottom_edge[:, 0], bottom_edge[:, 1], "s--", color=bottom_color,
            linewidth=1.8, markersize=5.0,
            label=f"Fault {index} bottom ({bot:g} km)", zorder=4,
        )
        top_handles.append(top_handle)
        bottom_handles.append(bottom_handle)

        ax_3d.add_collection3d(Poly3DCollection(
            [vertices_3d], facecolors=face_color, edgecolors="0.25",
            linewidths=1.0, alpha=0.55,
        ))
        top_vertices = vertices_3d[:2]
        bottom_vertices = vertices_3d[2:][::-1]
        ax_3d.plot(
            top_vertices[:, 0], top_vertices[:, 1], top_vertices[:, 2],
            "o-", color=top_color, linewidth=2.0, markersize=4.5,
            label=f"Fault {index} top ({top:g} km)",
        )
        ax_3d.plot(
            bottom_vertices[:, 0], bottom_vertices[:, 1], bottom_vertices[:, 2],
            "s--", color=bottom_color, linewidth=1.8, markersize=4.0,
            label=f"Fault {index} bottom ({bot:g} km)",
        )
        for top_point, bottom_point in zip(top_vertices, bottom_vertices):
            ax_3d.plot(
                [top_point[0], bottom_point[0]],
                [top_point[1], bottom_point[1]],
                [top_point[2], bottom_point[2]],
                color="0.35", linestyle="--", linewidth=0.8,
            )

    ax_map.set_title("Horizontal projection")
    ax_map.set_xlabel("East (km)")
    ax_map.set_ylabel("North (km)")
    ax_map.set_aspect("equal")
    ax_map.grid(linewidth=0.4, alpha=0.6)
    set_horizontal_limits(ax_map, projected)
    ax_map.legend(handles=top_handles + bottom_handles, loc="best", fontsize=8)

    depth_max = max(fault[-1] for fault in faults)
    minimum = projected.min(axis=0)
    maximum = projected.max(axis=0)
    span = np.maximum(maximum - minimum, 1.0)
    padding = 0.12 * max(span)
    ax_3d.set_title("3D finite-fault geometry")
    ax_3d.set_xlabel("East (km)", labelpad=8)
    ax_3d.set_ylabel("North (km)", labelpad=8)
    ax_3d.set_zlabel("Depth (km)", labelpad=8)
    ax_3d.set_xlim(minimum[0] - padding, maximum[0] + padding)
    ax_3d.set_ylim(minimum[1] - padding, maximum[1] + padding)
    ax_3d.set_zlim(0.0, max(depth_max * 1.12, 1.0))
    ax_3d.invert_zaxis()
    ax_3d.set_aspect("equal")
    ax_3d.view_init(elev=24.0, azim=-58.0)
    ax_3d.legend(loc="upper left", fontsize=8)

    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    plot_fault_geometry("fault.inp", "fault.svg")
