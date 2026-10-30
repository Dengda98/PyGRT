"""
    :file:     tensor.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2024-07-24

    动态和静态结果的张量计算、断层面投影及库仑应力

"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Sequence, Union

from .cli import _resolve_rcv_points, format_float, run_grt


PathLike = Union[str, os.PathLike]


def _run_static_file_module(path: PathLike, module: str, options: Sequence[object]) -> None:
    """运行只处理静态 NetCDF 文件的模块"""
    path = Path(path)
    if path.is_dir():
        raise ValueError(f"Only static synthesis files are supported: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Synthesis result does not exist: {path}")

    run_grt([module, f"-G{path}", *options])


def static_sproj(
    path: PathLike,
    *,
    strike: Optional[float] = None,
    dip: Optional[float] = None,
    rake: Optional[float] = None,
    rcv_points: Optional[PathLike] = None,
    force_rake: bool = False,
    **kwargs,
) -> None:
    """
    Project static stress tensors onto receiver-fault geometry in place.

    The input ``path`` must be a static synthesis NetCDF file, not a dynamic
    synthesis directory, and must contain the
    six stress components produced by ``static_stress``. For grid and ordinary
    points layouts, pass ``strike``, ``dip`` and ``rake`` together. For the
    finite receiver ``faults`` layout, pass only ``rake`` when rake values are undefined;
    set ``force_rake=True`` to replace every rake. ``rcv_points`` corresponds
    to the C module's ``-Q`` option and must contain six columns per row.

    Results are written back to ``path`` as ``sigma_n`` and ``tau_s``.

    :param    path:          Static synthesis NetCDF file containing stress components.
    :param    strike:        Manual receiver strike in degrees.
    :param    dip:           Manual receiver dip in degrees.
    :param    rake:          Manual receiver rake in degrees.
    :param    rcv_points:    Six-column receiver geometry file for the ``-Q`` option.
    :param    force_rake:    If true, append ``+f`` and force the manual rake for all finite points.
    """
    rcv_points = _resolve_rcv_points(rcv_points, kwargs, "static_sproj")
    options = []
    if strike is not None or dip is not None:
        if strike is None or dip is None or rake is None:
            raise ValueError("strike, dip and rake must be supplied together.")
        if force_rake:
            raise ValueError("force_rake requires rake alone.")
        geometry_text = "/".join(format_float(value) for value in (strike, dip, rake))
        options.append(f"-M{geometry_text}")
    elif rake is not None:
        geometry_text = format_float(rake)
        options.append(f"-M{geometry_text}{'+f' if force_rake else ''}")
    elif force_rake:
        raise ValueError("force_rake requires rake.")

    if rcv_points is not None:
        options.append(f"-Q{Path(rcv_points)}")

    _run_static_file_module(path, "static_sproj", options)


def sproj(path: PathLike, *, strike: Optional[float] = None, dip: Optional[float] = None,
          rake: Optional[float] = None, rcv_points: Optional[PathLike] = None,
          force_rake: bool = False) -> None:
    """
    Project dynamic stress tensors onto receiver-fault geometry.

    Accepts a single receiver directory or a multi-receiver root. Geometry is
    read from SAC receiver headers unless supplied explicitly. Finite receivers
    accept only ``rake``; ``force_rake`` overrides their existing rake values.
    With ``rcv_points``, the point count and coordinates must match in receiver-directory
    index order. All points are checked before any result is written.
    Writes ``sigma_n.sac`` and ``tau_s.sac`` in dyne/cm².

    :param    path:          Dynamic stress result directory.
    :param    strike:        Manual receiver strike in degrees.
    :param    dip:           Manual receiver dip in degrees.
    :param    rake:          Manual receiver rake in degrees.
    :param    rcv_points:    Six-column receiver geometry file in receiver-directory index order;
                             supply one point for a single receiver directory.
    :param    force_rake:    Force a manual rake for finite receivers.
    """
    path = Path(path)
    command = ["sproj", f"-G{path}"]
    # 显式机制和接收文件作为覆盖选项，缺省时由 C 模块读取 SAC 头段
    if strike is not None or dip is not None:
        if strike is None or dip is None or rake is None:
            raise ValueError("strike, dip and rake must be supplied together.")
        if force_rake:
            raise ValueError("force_rake requires rake alone.")
        command.append("-M" + "/".join(format_float(value) for value in (strike, dip, rake)))
    elif rake is not None:
        command.append(f"-M{format_float(rake)}{'+f' if force_rake else ''}")
    elif force_rake:
        raise ValueError("force_rake requires rake.")
    if rcv_points is not None:
        command.append(f"-Q{Path(rcv_points)}")
    run_grt(command)
    return None


def coulomb(path: PathLike, friction: float) -> None:
    """
    Compute dynamic Coulomb stress change, ``tau_s + friction * sigma_n``.

    Writes ``coulomb.sac`` in each receiver directory, in dyne/cm².

    :param    path:          Single receiver directory or multi-receiver root.
    :param    friction:      Nonnegative effective friction coefficient.
    """
    path = Path(path)
    run_grt(["coulomb", f"-G{path}", f"-F{format_float(friction)}"])
    return None


def compute_sproj(*args, **kwargs):
    """Legacy interface renamed to :func:`static_sproj`; calling it raises an error."""
    raise RuntimeError("compute_sproj() has been renamed to static_sproj(); use static_sproj() instead.")


def static_coulomb(path: PathLike, friction: float) -> None:
    """
    Compute Coulomb stress change in a static synthesis NetCDF file.

    The input file must already contain ``sigma_n`` and ``tau_s``, normally
    produced by :func:`static_sproj`. The result ``coulomb`` is written back
    to the same file using ``tau_s + friction * sigma_n``.

    :param    path:          Static synthesis NetCDF file containing ``sigma_n`` and ``tau_s``.
    :param    friction:      Nonnegative dimensionless effective friction coefficient.
    """
    _run_static_file_module(path, "static_coulomb", [f"-F{format_float(friction)}"])


def compute_coulomb(*args, **kwargs):
    """Legacy interface renamed to :func:`static_coulomb`; calling it raises an error."""
    raise RuntimeError("compute_coulomb() has been renamed to static_coulomb(); use static_coulomb() instead.")


def _run_dynamic_file_module(path: PathLike, module: str) -> None:
    """运行动态 SAC 目录后处理模块"""
    path = Path(path)
    if path.is_file():
        raise ValueError(f"Only dynamic synthesis directories are supported: {path}")
    if not path.is_dir():
        raise FileNotFoundError(f"Synthesis result does not exist: {path}")

    run_grt([module, f"-G{path}"])


def _run_static_tensor_module(path: PathLike, module: str) -> None:
    """运行静态 NetCDF 张量后处理模块"""
    path = Path(path)
    if path.is_dir():
        raise ValueError(f"Only static synthesis files are supported: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Synthesis result does not exist: {path}")

    run_grt([module, f"-G{path}"])


def strain(path: PathLike) -> None:
    """
    Compute a dynamic strain tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same SAC directory.

    :param    path:               Dynamic SAC synthesis directory.
    """
    _run_dynamic_file_module(path, "strain")


def static_strain(path: PathLike) -> None:
    """
    Compute a static strain tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same NetCDF file.

    :param    path:               Static synthesis NetCDF file.
    """
    _run_static_tensor_module(path, "static_strain")


def rotation(path: PathLike) -> None:
    """
    Compute a dynamic rotation tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same SAC directory.

    :param    path:               Dynamic SAC synthesis directory.
    """
    _run_dynamic_file_module(path, "rotation")


def static_rotation(path: PathLike) -> None:
    """
    Compute a static rotation tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same NetCDF file.

    :param    path:               Static synthesis NetCDF file.
    """
    _run_static_tensor_module(path, "static_rotation")


def stress(path: PathLike) -> None:
    """
    Compute a dynamic stress tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same SAC directory. Stress unit is dyne/cm² (= 0.1 Pa).

    :param    path:               Dynamic SAC synthesis directory.
    """
    _run_dynamic_file_module(path, "stress")


def static_stress(path: PathLike) -> None:
    """
    Compute a static stress tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same NetCDF file. Stress unit is dyne/cm² (= 0.1 Pa).

    :param    path:               Static synthesis NetCDF file.
    """
    _run_static_tensor_module(path, "static_stress")


def compute_strain(*args, **kwargs):
    """Legacy interface split into :func:`strain` and :func:`static_strain`; calling it raises an error."""
    raise RuntimeError("compute_strain() was replaced by strain() or static_strain(); use the matching interface instead.")


def compute_rotation(*args, **kwargs):
    """Legacy interface split into :func:`rotation` and :func:`static_rotation`; calling it raises an error."""
    raise RuntimeError("compute_rotation() was replaced by rotation() or static_rotation(); use the matching interface instead.")


def compute_stress(*args, **kwargs):
    """Legacy interface split into :func:`stress` and :func:`static_stress`; calling it raises an error."""
    raise RuntimeError("compute_stress() was replaced by stress() or static_stress(); use the matching interface instead.")
