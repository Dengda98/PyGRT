"""
    :file:     noise.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2026-09-18

    背景噪声模拟的点坐标生成和输入文件辅助函数

"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np

__all__ = ["random_points", "write_sources", "write_receivers"]


def random_points(
    count: int,
    *,
    north_range: Optional[Sequence[float]] = None,
    east_range: Optional[Sequence[float]] = None,
    radius_range: Optional[Sequence[float]] = None,
    center: Sequence[float] = (0.0, 0.0),
    depth: float = 0.0,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Generate local North-East-Depth points uniformly by area, in km.

    Choose a rectangle with ``north_range`` and ``east_range``, or an annulus with
    ``radius_range``. The returned array has shape ``(count, 3)``.

    :param    count:        Positive number of points.
    :param    north_range:  Rectangle north-coordinate limits, in km.
    :param    east_range:   Rectangle east-coordinate limits, in km.
    :param    radius_range: Inner and outer annulus radii, in km.
    :param    center:       Annulus center as north/east coordinates, in km.
    :param    depth:        Common nonnegative depth, in km.
    :param    seed:         Optional NumPy random seed for coordinate generation.
    """
    if isinstance(count, bool) or int(count) != count or count <= 0:
        raise ValueError("count must be a positive integer.")
    if not np.isfinite(depth) or depth < 0:
        raise ValueError("depth must be nonnegative and finite.")

    def limits(values, name):
        array = np.asarray(values, dtype=float)
        if array.shape != (2,) or not np.all(np.isfinite(array)) or array[0] > array[1]:
            raise ValueError(f"{name} must contain two finite ascending limits.")
        return array

    generator = np.random.default_rng(seed)
    result = np.full((int(count), 3), depth, dtype=float)
    if radius_range is not None:
        if north_range is not None or east_range is not None:
            raise ValueError("radius_range cannot be combined with rectangle limits.")
        radii = limits(radius_range, "radius_range")
        origin = np.asarray(center, dtype=float)
        if radii[0] < 0 or origin.shape != (2,) or not np.all(np.isfinite(origin)):
            raise ValueError("Annulus radii must be nonnegative and center must contain two finite values.")
        # 半径平方均匀，保证环带内按面积均匀分布
        radius = np.sqrt(generator.uniform(radii[0] ** 2, radii[1] ** 2, int(count)))
        phi = generator.uniform(0, 2 * np.pi, int(count))
        result[:, 0] = origin[0] + radius * np.cos(phi)
        result[:, 1] = origin[1] + radius * np.sin(phi)
    else:
        if north_range is None or east_range is None:
            raise ValueError("Specify radius_range or both north_range and east_range.")
        result[:, 0] = generator.uniform(*limits(north_range, "north_range"), int(count))
        result[:, 1] = generator.uniform(*limits(east_range, "east_range"), int(count))
    return result


def _point_rows(points, ids, prefix, fmt):
    """校验坐标及 ID，并生成四列点文件行"""
    coordinates = np.asarray(points, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or coordinates.shape[0] == 0:
        raise ValueError("points must be a nonempty (n, 3) array.")
    if not np.all(np.isfinite(coordinates)) or np.any(coordinates[:, 2] < 0):
        raise ValueError("Coordinates must be finite and depths nonnegative.")
    names = [f"{prefix}{i + 1:03d}" for i in range(len(coordinates))] if ids is None else [str(value) for value in ids]
    if len(names) != len(coordinates) or len(set(names)) != len(names):
        raise ValueError("ids must contain one unique name per point.")
    for name in names:
        if not name or len(name.encode("utf-8")) > 255 or any(c.isspace() or c in "/\\#" for c in name) or name in {".", ".."}:
            raise ValueError("IDs must be nonempty tokens without path separators or #.")
    return [[name, *(format(value, fmt) for value in point)] for name, point in zip(names, coordinates)]


def _write_rows(path, rows, header=None, parameters=None):
    """写出可选参数头，并按各列最大宽度对齐点坐标"""
    widths = [max(len(row[column]) for row in rows) for column in range(len(rows[0]))]
    if header:
        widths = [max(width, len(label)) for width, label in zip(widths, header)]
    target = Path(path)
    with target.open("w", encoding="utf-8") as stream:
        if parameters:
            width = max(len(name) for name, value in parameters)
            for name, value in parameters:
                stream.write(f"{name:<{width}}  {value}\n")
            stream.write("\n")
        if header:
            stream.write("  ".join(label.ljust(width) if column == 0 else label.rjust(width)
                                    for column, (label, width) in enumerate(zip(header, widths))) + "\n")
        for row in rows:
            stream.write("  ".join(value.ljust(width) if column == 0 else value.rjust(width)
                                    for column, (value, width) in enumerate(zip(row, widths))) + "\n")
    return target


def write_sources(
    path,
    points: Sequence[Sequence[float]],
    *,
    interval: float,
    force_z_range: Sequence[float] = (-1.0, 1.0),
    force_azimuth_range: Sequence[float] = (0.0, 360.0),
    strength_range: Sequence[float] = (1e15, 2e15),
    seed: int = 0,
    ids: Optional[Sequence[str]] = None,
    fmt: str = ".17g",
) -> Path:
    """Write an ordered sampling header and four-column source coordinates for ``grt synnoise``.

    The C module samples unit-force directions uniformly in solid angle within the
    specified downward-component and azimuth limits, then samples strengths uniformly.
    Range limits are written as separate whitespace-delimited values in the header.
    Row ``i`` fires once at ``i * interval`` seconds, starting at zero.

    :param    path:                Output source-file path.
    :param    points:              North-East-Depth coordinates in km, shape ``(n, 3)``.
    :param    interval:            Positive firing interval in seconds.
    :param    force_z_range:       Ascending unit-force downward-component limits within [-1, 1].
    :param    force_azimuth_range: Ascending clockwise-from-North angle limits within [0, 360], in degrees.
    :param    strength_range:      Ascending nonnegative force-strength limits, in dyne.
    :param    seed:                Nonnegative integer seed.
    :param    ids:                 Optional source IDs; defaults to SRC001, SRC002, etc.
    :param    fmt:                 Python format for coordinate columns; header values retain full precision.
    """
    rows = _point_rows(points, ids, "SRC", fmt)
    # 优先写出简短十进制，必要时保留可完整还原浮点值的精度
    def header_value(value):
        short = format(value, ".15g")
        return short if float(short) == value else repr(float(value))

    parameters = []
    for name, values, lower, upper in (("FORCE_Z", force_z_range, -1, 1),
                                       ("FORCE_AZIMUTH", force_azimuth_range, 0, 360),
                                       ("STRENGTH", strength_range, 0, np.inf)):
        limits = np.asarray(values, dtype=float)
        if (limits.shape != (2,) or not np.all(np.isfinite(limits)) or limits[0] > limits[1]
                or limits[0] < lower or limits[1] > upper):
            raise ValueError(f"{name} must contain two ascending limits within [{lower}, {upper}].")
        parameters.append((name, f"{header_value(limits[0])} {header_value(limits[1])}"))
    if not np.isfinite(interval) or interval <= 0:
        raise ValueError("interval must be positive and finite.")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or not 0 <= seed < 2**64:
        raise ValueError("seed must be an unsigned 64-bit integer.")
    parameters.extend([("INTERVAL", header_value(interval)), ("SEED", str(seed))])
    return _write_rows(path, rows, ["# ID", "north", "east", "depth"], parameters)


def write_receivers(
    path,
    points: Sequence[Sequence[float]],
    *,
    dt: float,
    length: float,
    ids: Optional[Sequence[str]] = None,
    fmt: str = ".17g",
) -> Path:
    """Write a recording header and receiver ID/North/East/Depth rows for ``grt synnoise``.

    :param    path:   Output receiver-file path.
    :param    points: North-East-Depth coordinates in km, shape ``(n, 3)``.
    :param    dt:     Positive sampling interval in seconds.
    :param    length: Positive record length in seconds; ``length/dt`` must be integral.
    :param    ids:    Optional receiver IDs; defaults to STA001, STA002, etc.
    :param    fmt:    Python format for coordinate columns; header values retain full precision.
    """
    rows = _point_rows(points, ids, "STA", fmt)
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be positive and finite.")
    if not np.isfinite(length) or length <= 0:
        raise ValueError("length must be positive and finite.")
    count = length / dt
    if not np.isfinite(count) or not 1 <= count <= 2147483647 or abs(count - round(count)) > 1e-6:
        raise ValueError("length/dt must be an integer supported by SAC v6 npts.")
    parameters = [("DT", repr(float(dt))), ("LENGTH", repr(float(length)))]
    return _write_rows(path, rows, ["# ID", "north", "east", "depth"], parameters)
