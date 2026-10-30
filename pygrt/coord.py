"""
    :file:     coord.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2024-07-24

    局部直角坐标与地理坐标之间的转换

"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Union

from .cli import format_float, run_grt


PathLike = Union[str, os.PathLike]


def _run_coordinate_transform(
    module: str,
    ingrid: Optional[PathLike],
    qfile: Optional[PathLike],
    outgrid: Optional[PathLike],
    lat0: Optional[float],
    lon0: Optional[float],
) -> None:
    """运行坐标转换模块，并校验其文件和参考点参数"""
    if (ingrid is None) == (qfile is None):
        raise ValueError("Specify exactly one of ingrid and qfile.")
    if outgrid is None:
        raise ValueError("outgrid is required.")
    if lat0 is None or lon0 is None:
        raise ValueError("lat0 and lon0 are required.")

    try:
        lat0 = float(lat0)
        lon0 = float(lon0)
    except (TypeError, ValueError):
        raise ValueError("lat0 and lon0 must be numbers.") from None
    if not (-90.0 < lat0 < 90.0):
        raise ValueError("lat0 must be in (-90, 90).")
    if not (-180.0 <= lon0 <= 180.0):
        raise ValueError("lon0 must be in [-180, 180].")

    input_path = Path(ingrid if ingrid is not None else qfile)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input coordinate file does not exist: {input_path}")
    output_path = Path(outgrid)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    input_option = f"-G{input_path}" if ingrid is not None else f"-Q{input_path}"
    run_grt([
        module,
        input_option,
        f"-O{output_path}",
        f"-C{format_float(lat0)}/{format_float(lon0)}",
    ])


def xy2geo(
    ingrid: Optional[PathLike] = None,
    *,
    qfile: Optional[PathLike] = None,
    outgrid: Optional[PathLike] = None,
    lat0: Optional[float] = None,
    lon0: Optional[float] = None,
) -> None:
    """
    Convert local north/east coordinates to geographic latitude/longitude.

    Specify exactly one of ``ingrid`` and ``qfile``. ``ingrid`` is a static
    NetCDF input file, while ``qfile`` is a text coordinate file. For a text
    input, only the first two columns are converted and the remaining text is
    preserved. The reference point is passed as ``-Clat0/lon0``.

    :param    ingrid: Static NetCDF input file for the ``-G`` option.
    :param    qfile: Text coordinate input file for the ``-Q`` option.
    :param    outgrid: Output NetCDF or text file for the ``-O`` option.
    :param    lat0: Reference latitude in degree.
    :param    lon0: Reference longitude in degree.
    """
    _run_coordinate_transform("xy2geo", ingrid, qfile, outgrid, lat0, lon0)


def geo2xy(
    ingrid: Optional[PathLike] = None,
    *,
    qfile: Optional[PathLike] = None,
    outgrid: Optional[PathLike] = None,
    lat0: Optional[float] = None,
    lon0: Optional[float] = None,
) -> None:
    """
    Convert geographic latitude/longitude to local north/east coordinates.

    Specify exactly one of ``ingrid`` and ``qfile``. ``ingrid`` is a NetCDF
    input file with ``lat``/``lon`` coordinates, while ``qfile`` is a text
    coordinate file. For a text input, only the first two columns are
    converted and the remaining text is preserved. The reference point is
    passed as ``-Clat0/lon0``.

    :param    ingrid: NetCDF input file for the ``-G`` option.
    :param    qfile: Text coordinate input file for the ``-Q`` option.
    :param    outgrid: Output NetCDF or text file for the ``-O`` option.
    :param    lat0: Reference latitude in degree.
    :param    lon0: Reference longitude in degree.
    """
    _run_coordinate_transform("geo2xy", ingrid, qfile, outgrid, lat0, lon0)
