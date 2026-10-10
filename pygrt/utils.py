"""
    :file:     utils.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包含一些数据处理操作上的补充以及其他辅助函数

"""

from __future__ import annotations

import os
import glob
import warnings
from contextlib import contextmanager
from enum import IntFlag
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from scipy.interpolate import interpn
from scipy.io import netcdf_file
from scipy.special import jv
import numpy.ctypeslib as npct

from .cli import format_float, format_range, run_grt
from .c_interfaces import (
    C_grt_lamb_parse_phase_list,
    C_grt_solve_lamb1,
    C_grt_solve_lamb2,
    C_grt_solve_lamb3,
)


__all__ = [
    "sproj",
    "coulomb",
    "read_nc",
    "read_nc_variables",
    "okada",
    "lamb",
    "strain",
    "rotation",
    "stress",
    "static_strain",
    "static_rotation",
    "static_stress",
    "static_sproj",
    "static_coulomb",
    "compute_okada",
    "compute_strain",
    "compute_rotation",
    "compute_stress",
    "compute_sproj",
    "compute_coulomb",
    "xy2geo",
    "geo2xy",
    "read_kernels_freqs",
    "read_statsfile",
    "read_statsfile_ptam",
    "plot_statsdata",
    "plot_statsdata_ptam",
    "read_dispersion",
    "read_eigenfunction",
    "read_energy_integral",
    "read_sensitivity",
    "read_secfunc",
    "plot_dispersion",
    "plot_eigenfunction",
    "plot_sensitivity",
    "plot_secfunc",
    "lamb1",
    "solve_lamb1",
    "lamb2",
    "lamb3",
]


PathLike = Union[str, os.PathLike]


@contextmanager
def _temporary_array_option(values: Optional[np.ndarray], flag: str):
    """将多值数组通过临时文件传给 CLI 选项，退出上下文时自动清理"""
    if values is None:
        yield None
    elif values.size == 1:
        yield f"-{flag}{format_float(values[0])}"
    else:
        with TemporaryDirectory(prefix=f"pygrt_{flag}_") as tmpdir:
            path = Path(tmpdir) / "values.txt"
            # 每行一个数值，保持与直接传值时相同的数值精度
            np.savetxt(path, values, fmt="%.15g")
            yield f"-{flag}{path}"


def _resolve_rcv_points(rcv_points: Optional[PathLike], kwargs: dict, function_name: str):
    """解析 ``recv_points`` 兼容关键字并返回规范参数"""
    if "recv_points" in kwargs:
        legacy_points = kwargs.pop("recv_points")
        if rcv_points is not None and legacy_points is not None:
            raise TypeError(f"{function_name}() got both 'rcv_points' and deprecated 'recv_points'.")
        warnings.warn(
            "'recv_points' is deprecated; it is an alias of 'rcv_points'.",
            DeprecationWarning,
            stacklevel=3,
        )
        if rcv_points is None:
            rcv_points = legacy_points

    if kwargs:
        name = next(iter(kwargs))
        raise TypeError(f"{function_name}() got an unexpected keyword argument {name!r}.")
    return rcv_points


QWV_NUM = 3
INTEG_NUM = 4
SRC_M_NUM = 6
SRC_M_ORDERS = [0, 0, 1, 0, 1, 2]
SRC_M_NAME_ABBR = ["EX", "VF", "HF", "DD", "DS", "SS"]
qwvchs = ["q", "w", "v"]
NPCT_REAL_TYPE = "f8"
NPCT_CMPLX_TYPE = "c16"


def read_nc(path: PathLike) -> dict:
    """
    Read a NetCDF file into a nested dictionary.

    The returned dictionary contains three top-level entries:

    * ``dimensions`` - mapping from dimension name to length
    * ``variables`` - mapping from variable name to
      ``{"dimensions", "data", "attributes"}``
    * ``attributes`` - global NetCDF attributes

    Variable arrays are available at ``variables[name]["data"]``.

    :param    path:               Path to the NetCDF file.

    :return: A dictionary containing the NetCDF data and metadata.
    """
    def attribute_value(value):
        # 将 NetCDF 属性转为更方便使用的 Python 值
        if isinstance(value, np.ndarray) and value.ndim == 0:
            return value.item()
        if isinstance(value, bytes):
            return value.decode("utf-8")
        return value

    path = str(path)
    if not Path(path).is_file():
        raise FileNotFoundError(f"NetCDF file does not exist: {path}")

    with netcdf_file(path, mode="r", mmap=False) as dataset:
        dimensions = {name: int(length) for name, length in dataset.dimensions.items()}
        attributes = {name: attribute_value(getattr(dataset, name)) for name in dataset._attributes}
        variables = {}
        result = {
            "dimensions": dimensions,
            "variables": variables,
            "attributes": attributes,
        }
        for name, variable in dataset.variables.items():
            data = np.array(variable[:], copy=True)
            variable_attributes = {key: attribute_value(value) for key, value in variable._attributes.items()}
            variables[name] = {
                "dimensions": tuple(variable.dimensions),
                "data": data,
                "attributes": variable_attributes,
            }
    return result


def read_nc_variables(path: PathLike) -> dict:
    """
    Read NetCDF variables as a name-to-array mapping.

    This is a simplified form of :func:`read_nc`. Global attributes and
    dimension metadata are omitted.

    :param    path:               Path to the NetCDF file.

    :return: A dictionary mapping each variable name to a ``numpy.ndarray``.
    """
    return {name: info["data"] for name, info in read_nc(path)["variables"].items()}





def okada(
    *,
    modelparams: Sequence[float],
    depsrc: Optional[float] = None,
    deprcv: Optional[float] = None,
    norths: Optional[Sequence[float]] = None,
    easts: Optional[Sequence[float]] = None,
    rcv_points: Optional[PathLike] = None,
    rcv_fault: Optional[PathLike] = None,
    rcv_fault_size: Optional[Sequence[float]] = None,
    output_path: PathLike,
    scale: Optional[float] = None,
    scale_with_mu: bool = False,
    strike: Optional[float] = None,
    dip: Optional[float] = None,
    rake: Optional[float] = None,
    src_fault: Optional[PathLike] = None,
    zne: bool = False,
    calc_upar: bool = False,
    **kwargs,
) -> None:
    r"""
    Synthesize static displacement with the Okada homogeneous half-space solution.

    Results are written to the NetCDF file ``output_path``. The source, receiver
    grid, component and derivative arguments are intentionally close to
    :meth:`PyModel1D.static_syn`, but Okada only needs the homogeneous
    half-space model parameters ``(vp, vs, rho)`` and does not require a static
    Green's function file.

    The point-source type is inferred from the source-specific parameters:

    * no ``strike``, ``dip`` or ``rake`` - explosion (``EX``)
    * ``strike`` and ``dip`` - tensile crack (``TS``), or double-couple (``DC``)
      when ``rake`` is also supplied

    A Coulomb-format finite-fault file can be passed through ``src_fault``.
    Its Kode column selects the rectangular or point-source interpretation of
    the two slip columns. If the seventh header column is exactly ``rake``,
    the two values are interpreted as rake/net slip; the filename suffix is
    not used to select this format. The rake/net-slip interpretation supports
    Kode 100 only.
    The finite fault is evaluated directly as Okada rectangular patches, so
    no source subdivision option is needed.

    ``strike``, ``dip`` and ``rake`` must be supplied as a complete geometry
    when they are used. They are mutually exclusive with ``src_fault``.

    All arguments must be passed by keyword.

    :param    modelparams:      Homogeneous half-space parameters ``(vp, vs, rho)``;
                                velocities are in km/s and density is in g/cm^3
    :param    depsrc:           Point-source depth in km. Required for point
                                sources and forbidden for finite faults
    :param    deprcv:           Receiver depth in km for a regular grid. Forbidden
                                when ``rcv_points`` is used
    :param    norths:           North grid range ``(start, stop, step)`` in km
    :param    easts:            East grid range ``(start, stop, step)`` in km
    :param    rcv_points:       ASCII receiver file with either ``north east depth``
                                or ``north east depth strike dip rake``; coordinates
                                are in km and angles are in degrees
    :param    rcv_fault:        Coulomb-format finite receiver-fault file with Kode=100. Each
                                fault contributes its center, or subfault centers
                                when ``rcv_fault_size`` is supplied. Slip magnitude is ignored;
                                an exact ``rake`` header preserves the angle even at zero slip,
                                otherwise the slip columns define direction
    :param    rcv_fault_size:   Optional ``(dL, dW)`` receiver subdivision size
                                in km along strike / dip. Both values must be positive or zero.
                                Omit or set ``(0, 0)`` to use each fault center
    :param    output_path:      Output NetCDF file path
    :param    scale:            Point-source scale in dyne-cm unless
                                ``scale_with_mu`` is true. Not used for finite faults
    :param    scale_with_mu:    If true, pass ``-Su`` and treat ``scale`` as potency
                                or area times slip in cm^3
    :param    strike:           Fault strike in degrees, in [0, 360]
    :param    dip:              Fault dip in degrees, in [0, 90]
    :param    rake:             Slip rake in degrees, in [-180, 180]
    :param    src_fault:        Coulomb-format finite-fault file with 11 data columns;
                                an exact ``rake`` token in the seventh header column
                                selects Kode 100 rake/net-slip interpretation. Mutually
                                exclusive with point-source options
    :param    zne:              If true, output ZNE instead of ZRT components
    :param    calc_upar:        If true, also output spatial displacement derivatives
    """

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(modelparams, (str, bytes)):
        raise TypeError("modelparams must be a sequence of (vp, vs, rho).")
    try:
        if len(modelparams) != 3:
            raise ValueError("modelparams must contain exactly three values: (vp, vs, rho).")
    except TypeError:
        raise TypeError("modelparams must be a sequence of (vp, vs, rho).") from None
    rcv_points = _resolve_rcv_points(rcv_points, kwargs, "okada")
    vp, vs, rho = modelparams

    use_ff = src_fault is not None
    use_q = rcv_points is not None
    use_u = rcv_fault is not None
    use_xy = norths is not None or easts is not None
    has_strike = strike is not None
    has_dip = dip is not None
    has_rake = rake is not None
    has_mechanism = has_strike or has_dip or has_rake
    has_point_source_options = scale is not None or depsrc is not None or has_mechanism or scale_with_mu

    def source_option() -> Optional[str]:
        if not has_mechanism:
            return None
        if not has_strike or not has_dip:
            raise ValueError("strike and dip must be supplied together.")
        if has_rake:
            return f"-M{format_float(strike)}/{format_float(dip)}/{format_float(rake)}"
        return f"-M{format_float(strike)}/{format_float(dip)}"

    if ((use_q or use_u) and use_xy):
        raise ValueError("rcv_points/rcv_fault is mutually exclusive with norths/easts.")
    if (use_q and use_u):
        raise ValueError("rcv_points and rcv_fault are mutually exclusive.")
    if use_xy and (norths is None or easts is None):
        raise ValueError("norths and easts must be supplied together.")
    if ((not use_q) and (not use_u) and (not use_xy)):
        raise ValueError("Specify norths/easts, rcv_points or rcv_fault.")
    if depsrc is not None and depsrc < 0.0:
        raise ValueError("depsrc must be nonnegative.")
    if deprcv is not None and deprcv < 0.0:
        raise ValueError("deprcv must be nonnegative.")
    if ((use_q or use_u) and (deprcv is not None)):
        raise ValueError("rcv_points/rcv_fault is mutually exclusive with deprcv.")
    if rcv_fault_size is not None:
        if (not use_u):
            raise ValueError("rcv_fault_size requires rcv_fault.")
        if len(rcv_fault_size) != 2 or not (all(value == 0 for value in rcv_fault_size) or all(value > 0 for value in rcv_fault_size)):
            raise ValueError("rcv_fault_size must contain two zeros or two positive values.")

    command = [
        "okada",
        f"-H{format_float(vp)}/{format_float(vs)}/{format_float(rho)}",
        f"-O{output}",
    ]
    if use_ff:
        if has_point_source_options:
            raise ValueError("src_fault is mutually exclusive with point-source options.")
        command.append(f"-C{Path(src_fault)}")
    else:
        if scale is None:
            raise ValueError("scale is required for point-source synthesis.")
        if depsrc is None:
            raise ValueError("depsrc is required for point-source synthesis.")
        command.append(f"-S{'u' if scale_with_mu else ''}{format_float(scale)}")
        command.append(f"-Ds{format_float(depsrc)}")
        geometry_option = source_option()
        if geometry_option is not None:
            command.append(geometry_option)

    if deprcv is not None:
        command.append(f"-Dr{format_float(deprcv)}")
    elif (not use_q) and (not use_u):
        raise ValueError("deprcv is required for grid receivers.")

    if use_q:
        command.append(f"-Q{Path(rcv_points)}")
    elif use_u:
        receiver_option = f"-U{Path(rcv_fault)}"
        if rcv_fault_size is not None:
            receiver_option += f"+i{format_float(rcv_fault_size[0])}/{format_float(rcv_fault_size[1])}"
        command.append(receiver_option)
    else:
        command.append(f"-X{format_range(norths, 'norths')}")
        command.append(f"-Y{format_range(easts, 'easts')}")

    if zne:
        command.append("-N")
    if calc_upar:
        command.append("-e")

    run_grt(command)
    return None


def compute_okada(*args, **kwargs):
    """Legacy interface renamed to :func:`okada`; calling it raises an error."""
    raise RuntimeError("compute_okada() has been renamed to okada(); use okada() instead.")


def _run_static_file_module(path: PathLike, module: str, options: Sequence[object]) -> None:
    """运行只处理静态 NetCDF 文件的模块"""
    path = Path(path)
    if path.is_dir():
        raise ValueError(f"Only static synthesis files are supported: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Synthesis result does not exist: {path}")

    run_grt([module, f"-G{path}", *options])


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

    run_grt([module, path])


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


#=================================================================================================================
#
#                                           积分过程文件读取及绘制
#
#=================================================================================================================


def read_statsfile(statsfile:str):
    '''
        read a statsfile  

        :param    statsfile:       File path (Wildcards can be used to simplify input)

        :return:
            - **data** -     `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ custom type array 
    '''
    Lst = glob.glob(statsfile)
    if len(Lst) != 1:
        raise OSError(f"{statsfile} should only match one file, but {len(Lst)} matched.")
    statsfile = Lst[0]
    print(f"read in {statsfile}.")

    basename = os.path.basename(statsfile)

    # 确定自定义数据类型  EX_q, EX_w, VF_q, ...
    dtype = [('k' if basename[0] == 'K' else 'c', NPCT_REAL_TYPE)]
    for im in range(SRC_M_NUM):
        modr = SRC_M_ORDERS[im]
        for c in range(QWV_NUM):
            if modr==0 and qwvchs[c] == 'v':
                continue 

            dtype.append((f"{SRC_M_NAME_ABBR[im]}_{qwvchs[c]}", NPCT_CMPLX_TYPE))


    data = np.fromfile(statsfile, dtype=dtype)

    return data


def read_kernels_freqs(statsdir:str, vels:Union[np.ndarray,None]=None, ktypes:Union[List[str],None]=None):
    r"""
        read all statsfiles in statsdir (except that of 0 frequency).
        If record wavenumber, interpolate to the phase velocity.

        :param        statsdir:     directory path
        :param        vels:         When a positive-order vels (km/s) is specified, files starting with `K_` are read 
                                    and linear interpolation from wavenumber to phase velocity is performed.
                                    Otherwise read the files starting with `C_`
        :param        ktypes:       Specify the return of a series of kernel function names,
                                    such as `EX_q`, `DS_w`, etc. By default, all are returned

        :return:
            - **kerDct**  -   kernel functions in a dict
    """

    dointerp = vels is not None

    if (dointerp) and not np.all(np.diff(vels) > 0):
        raise ValueError("vels must be in ascending order.")
    
    K_statspaths = glob.glob(os.path.join(statsdir, "K_*"))
    if len(K_statspaths) == 0 and dointerp:
        raise ValueError("You want to interpolate from k to c, but found 0 statsfiles recording k.")
    
    C_statspaths = glob.glob(os.path.join(statsdir, "C_*"))
    if len(C_statspaths) == 0 and not dointerp:
        raise ValueError("Found 0 statsfiles directly recording c.")
    
    statspaths = K_statspaths if dointerp else C_statspaths

    KLst = np.array(statspaths)
    freqs = np.array([float(s.split("_")[-1]) for s in KLst])
    # 根据freqs排序
    _idx = np.argsort(freqs)
    freqs[:] = freqs[_idx]
    KLst[:] = KLst[_idx]
    del _idx 

    # 去除零频
    if freqs[0] == 0.0:
        freqs = freqs[1:]
        KLst = KLst[1:]

    kerDct = {}
    kerDct['_vels'] = vels.copy() if dointerp else []
    kerDct['_freqs'] = freqs.copy()

    for i in range(len(freqs)):
        Kpath = KLst[i]
        freq = freqs[i]
        w = 2*np.pi*freq

        data = read_statsfile(Kpath)
        
        if dointerp:
            v = w/data['k']

            # 检查v范围
            v1 = np.min(v)
            v2 = np.max(v)
            if v1 > vels[0] or v2 < vels[-1]:
                raise ValueError(f"In freq={freq:.5e}, minV={v1:.5e}, maxV={v2:.5e}, insufficient wavenumber samples"
                                " to interpolate on vels.")
        else:
            if len(kerDct['_vels']) == 0:
                kerDct['_vels'] = data['c'].copy()

        for key in data.dtype.names:
            if key == 'k' or key == 'c':
                continue 
            if (ktypes is not None) and (key not in ktypes):
                continue 

            if key not in kerDct.keys():
                kerDct[key] = []

            if dointerp:
                # 如果越界会报错
                F = interpn((v,), data[key], vels)
                kerDct[key].append(F)
            else:
                kerDct[key].append(data[key])

    # 将每个核函数结果拼成2D数组
    for key in kerDct.keys():
        if key[0] == '_':
            continue
        kerDct[key] = np.vstack(kerDct[key])

    return kerDct


def read_statsfile_ptam(statsfile:str):
    '''
        read a statsfile from PTAM process  

        :param    statsfile:       PTAM stats file path. The basename must start
                                   with ``PTAM``. Wildcards can be used to
                                   simplify input.

        :return:
            - **data1** -     `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ custom type array, during DCM or (SA)FIM
            - **data2** -     `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ custom type array, during PTAM
            - **ptam_data** -   `numpy.ndarray <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>`_ custom type array, record the peak/trough from PTAM
            - **dist** -      epicentral distance from the filename (km)
    '''
    Lst = glob.glob(statsfile)
    if len(Lst) != 1:
        raise OSError(f"{statsfile} should only match one file, but {len(Lst)} matched.")
    statsfile = Lst[0]

    # 必须是 PTAM 开头的统计文件，避免误读 K/C 核函数文件
    PTAMname = os.path.basename(statsfile)
    if not PTAMname.startswith("PTAM"):
        raise ValueError(
            f"{statsfile} is not a PTAM stats file; the basename must start with 'PTAM'."
        )

    # 获得震中距
    dist = float(os.path.dirname(statsfile).split("_")[-1])

    # 从文件路径命名中，获得对应的K文件路径
    if "_" in PTAMname:  # 动态解
        splits = PTAMname.split("_")
        splits[-3] = "K"
        K_basename= "_".join(splits)
    else:
        K_basename = "K" # 静态解
        
    data1 = read_statsfile(os.path.join(os.path.dirname(os.path.dirname(statsfile)), K_basename))
    data2 = read_statsfile(os.path.join(os.path.dirname(statsfile), K_basename))

    # 确定自定义数据类型  sum_EX_0_k, sum_EX_0, sum_VF_0_k, ...
    # 各格林函数数值积分的值(k上限位于不同的波峰波谷)
    # 开头的sum表示这是波峰波谷位置处的数值积分的值(不含dk)，
    # 末尾的k表示对应积分值的波峰波谷位置的k值
    dtype = []
    for im in range(SRC_M_NUM):
        modr = SRC_M_ORDERS[im]
        for v in range(INTEG_NUM):
            if modr==0 and v!=0 and v!=2:
                continue 

            dtype.append((f"sum_{SRC_M_NAME_ABBR[im]}_{v}_k", NPCT_REAL_TYPE))
            dtype.append((f"sum_{SRC_M_NAME_ABBR[im]}_{v}", NPCT_CMPLX_TYPE))


    ptam_data = np.fromfile(statsfile, dtype=dtype)

    return data1, data2, ptam_data, dist



def _get_stats_Fname(statsdata:np.ndarray, karr:np.ndarray, dist:float, srctype:str, ptype:str):
    # 根据ptype获得对应的核函数
    krarr = karr*dist

    # 从数组中找到震源名称的索引
    try:
        _idx = SRC_M_NAME_ABBR.index(srctype)
        mtype = str(SRC_M_ORDERS[_idx])
    except:
        raise ValueError(f"{srctype} is an invalid name.")

    if mtype=='0':
        if ptype=='0':
            Fname = rf"$F(k,\omega)=q^{{({srctype})}}(k, \omega)$"
            Farr = statsdata[f'{srctype}_q']
            FJname = rf"$ - F(k,\omega)J_1(kr)k$"
            FJarr =  - jv(1, krarr) * Farr * karr
        elif ptype=='2':
            Fname = rf"$F(k,\omega)=w^{{({srctype})}}(k, \omega)$"
            FJname = rf"$F(k,\omega)J_0(kr)k$"
            Farr = statsdata[f'{srctype}_w']
            FJarr = jv(0, krarr) * Farr * karr
        else:
            raise ValueError(f"source {srctype}, m={mtype}, p={ptype} is not supported.")
        
    elif mtype in ['1', '2']:
        m = int(mtype)
        if ptype=='0':
            Fname = rf"$F(k,\omega)=q^{{({srctype})}}(k, \omega)$"
            Farr = statsdata[f'{srctype}_q']
            FJname = rf"$F(k,\omega)J_{m-1}(kr)k$"
            FJarr = jv(m-1, krarr) * Farr * karr
        elif ptype=='1':
            Fname = rf"$F(k,\omega)=q^{{({srctype})}}(k, \omega) + v^{{({srctype})}}(k, \omega)$"
            Farr = (statsdata[f'{srctype}_q'] + statsdata[f'{srctype}_v'])
            FJname = rf"$ - F(k,\omega) \dfrac{{{m}}}{{kr}} J_{m}(kr)k$"
            FJarr =  - jv(m, krarr) * Farr * m/dist
        elif ptype=='2':
            Fname = rf"$F(k,\omega)=w^{{({srctype})}}(k, \omega)$"
            Farr = statsdata[f'{srctype}_w']
            FJname = rf"$F(k,\omega)J_{m}(kr)k$"
            FJarr = jv(m, krarr) * Farr * karr
        elif ptype=='3':
            Fname = rf"$F(k,\omega)=v^{{({srctype})}}(k, \omega)$"
            Farr = statsdata[f'{srctype}_v']
            FJname = rf"$ - F(k,\omega)J_{m-1}(kr)k$"
            FJarr =  - jv(m-1, krarr) * Farr * karr
        else:
            raise ValueError(f"source {srctype}, m={mtype}, p={ptype} is not supported.")
        
    else:
        raise ValueError(f"source {srctype}, m={mtype}, p={ptype} is not supported.")
    
    return Fname, Farr, FJname, FJarr


def plot_statsdata(statsdata:np.ndarray, dist:float, srctype:str, ptype:str, RorI:Union[bool,int]=True,
                   fig:Union[Figure,None]=None, axs:Union[Axes,None]=None):
    r'''
        Based on the data read by the :func:`read_statsfile <pygrt.utils.read_statsfile>` function,
        plot the kernel function :math:`F(k,\omega)`, the integrand :math:`F(k,\omega)J_m(kr)k`, 
        and calculate the cumulative integral :math:`\sum F(k,\omega)J_m(kr)k` .

        .. note:: Not every source type corresponds to every order and every integration type, see :ref:`grn_types` for details.

        :param    statsdata:         return value of :func:`read_statsfile <pygrt.utils.read_statsfile>` function
        :param    dist:              epicentral distance (km)
        :param    srctype:           abbreviation of source type, including EX, VF, HF, DD, DS, SS
        :param    ptype:             integration type (0,1,2,3)
        :param    RorI:              whether to plot real or imaginary part, default is real part, pass 2 to plot both
        :param    fig:               user-defined matplotlib.Figure object, default is None
        :param    axs:               user-defined matplotlib.Axes object array (three elements), default is None

        :return:
                - **fig** -                        matplotlib.Figure object
                - **(ax1,ax2,ax3)** -              matplotlib.Axes object array
    '''

    ptype = str(ptype)

    karr = statsdata['k'] 
    dk = (karr[1] - karr[0])   # 假设均匀dk
    is_evendk = np.allclose(np.diff(karr), dk, atol=1e-10)  # 是否为均匀dk
    if not is_evendk:
        raise ValueError("Sorry, this function only supports even-distributed k.")
    
    if 0.5*np.pi/dk < dist:  # 对于bessel函数这种震荡函数，假设一个周期内至少取4个点
        print(f"WARNING! dist ({dist}) > PI/(2*dk) ({0.5*np.pi/dk:.5e}.)")

    Fname, Farr, FJname, FJarr = _get_stats_Fname(statsdata, karr, dist, srctype, ptype)
    
    if fig is None or axs is None:
        fig, axs = plt.subplots(3, 1, figsize=(8, 9), gridspec_kw=dict(hspace=0.7))
    
    # axs长度必须为三个
    if len(axs) != 3:
        raise ValueError("axs should have 3 elements.")

    ax1, ax2, ax3 = axs

    if isinstance(RorI, int) and RorI==2:
        ax1.plot(karr, np.real(Farr), lw=0.8, label='Real') 
        ax1.plot(karr, np.imag(Farr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax1.plot(karr, np.real(Farr), lw=0.8, label='Real') 
        else:
            ax1.plot(karr, np.imag(Farr), lw=0.8, label='Imag') 

    ax1.set_xlabel('k /$km^{-1}$')
    ax1.set_title(Fname)
    ax1.grid()
    ax1.legend(loc='lower left')

    if isinstance(RorI, int) and RorI==2:
        ax2.plot(karr, np.real(FJarr), lw=0.8, label='Real') 
        ax2.plot(karr, np.imag(FJarr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax2.plot(karr, np.real(FJarr), lw=0.8, label='Real') 
        else:
            ax2.plot(karr, np.imag(FJarr), lw=0.8, label='Imag') 
    ax2.set_title(FJname)
    ax2.set_xlabel('k /$km^{-1}$')
    ax2.grid()
    ax2.legend(loc='lower left')

    # 数值积分，不乘系数dk 
    Parr = np.cumsum(FJarr)

    if isinstance(RorI, int) and RorI==2:
        ax3.plot(karr, np.real(Parr), lw=0.8, label='Real') 
        ax3.plot(karr, np.imag(Parr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax3.plot(karr, np.real(Parr), lw=0.8, label='Real') 
        else:
            ax3.plot(karr, np.imag(Parr), lw=0.8, label='Imag') 
    ax3.set_title(rf'$\sum_k$ {FJname}')
    ax3.set_xlabel("k /$km^{-1}$")
    ax3.grid()
    ax3.legend(loc='lower left')

    return fig, (ax1, ax2, ax3)


def plot_statsdata_ptam(statsdata1:np.ndarray, statsdata2:np.ndarray, statsdata_ptam:np.ndarray,
                        dist:float, srctype:str, ptype:str, RorI:Union[bool,int]=True,
                        fig:Union[Figure,None]=None, axs:Union[Axes,None]=None):
    r'''
        Based on data read by the :func:`read_statsfile_ptam <pygrt.utils.read_statsfile_ptam>` function,
        simply calculate and plot the cumulative integral as well as the peak/trough positions used by PTAM.

        .. note:: Not every source type corresponds to every order and every integration type, see :ref:`grn_types` for details.

        :param    statsdata1:        integral process data during DWM or FIM
        :param    statsdata2:        integral process data during PTAM
        :param    statsdata_ptam:    peak/trough positions and amplitudes from PTAM
        :param    dist:              epicentral distance (km)
        :param    srctype:           abbreviation of source type, including EX, VF, HF, DD, DS, SS  
        :param    ptype:             integration type (0, 1, 2, 3)
        :param    RorI:              whether to plot real or imaginary part, default is real part, pass 2 to plot both
        :param    fig:               user-defined matplotlib.Figure object, default is None
        :param    axs:               user-defined matplotlib.Axes object array (three elements), default is None

        :return:  
                - **fig** -                        matplotlib.Figure object   
                - **(ax1, ax2, ax3)** -            matplotlib.Axes object array
    '''

    ptype = str(ptype)

    karr1 = statsdata1['k'] 
    dk1 = karr1[1] - karr1[0]
    Fname, Farr1, FJname, FJarr1 = _get_stats_Fname(statsdata1, karr1, dist, srctype, ptype)
    karr2 = statsdata2['k'] 
    dk2 = karr2[1] - karr2[0]
    Fname, Farr2, FJname, FJarr2 = _get_stats_Fname(statsdata2, karr2, dist, srctype, ptype)

    is_evendk = np.allclose(np.diff(karr1), dk1, atol=1e-10) and np.allclose(np.diff(karr2), dk2, atol=1e-10)  # 是否为均匀dk
    if not is_evendk:
        raise ValueError("Sorry, this function only supports even-distributed k.")

    # 将两个过程的结果拼起来
    Farr = np.hstack((Farr1, Farr2))
    karr = np.hstack((karr1, karr2))
    FJarr = np.hstack((FJarr1, FJarr2))

    if fig is None or axs is None:
        fig, axs = plt.subplots(3, 1, figsize=(8, 9), gridspec_kw=dict(hspace=0.7))
    
    # axs长度必须为三个
    if len(axs) != 3:
        raise ValueError("axs should have 3 elements.")

    ax1, ax2, ax3 = axs

    if isinstance(RorI, int) and RorI==2:
        ax1.plot(karr, np.real(Farr), lw=0.8, label='Real') 
        ax1.plot(karr, np.imag(Farr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax1.plot(karr, np.real(Farr), lw=0.8, label='Real') 
        else:
            ax1.plot(karr, np.imag(Farr), lw=0.8, label='Imag') 

    ax1.set_xlabel('k /$km^{-1}$')
    ax1.set_title(Fname)
    ax1.grid()
    ax1.legend(loc='lower left')

    if isinstance(RorI, int) and RorI==2:
        ax2.plot(karr, np.real(FJarr), lw=0.8, label='Real') 
        ax2.plot(karr, np.imag(FJarr), lw=0.8, label='Imag') 
    else:
        if RorI:
            ax2.plot(karr, np.real(FJarr), lw=0.8, label='Real') 
        else:
            ax2.plot(karr, np.imag(FJarr), lw=0.8, label='Imag') 
    ax2.set_title(FJname)
    ax2.set_xlabel('k /$km^{-1}$')
    ax2.grid()
    ax2.legend(loc='lower left')

    # 波峰波谷位置，用红十字标记
    ptKarr = statsdata_ptam[f'sum_{srctype}_{ptype}_k']
    ptFJarr = statsdata_ptam[f'sum_{srctype}_{ptype}']

    # 数值积分，不乘系数dk 
    Parr1 = np.cumsum(FJarr1) 
    Parr2 = np.cumsum(FJarr2)  
    Parr = np.hstack([Parr1, Parr2*dk2/dk1+Parr1[-1]])

    if isinstance(RorI, int) and RorI==2:
        ax3.plot(karr, np.real(Parr), lw=0.8, label='Real') 
        ax3.plot(ptKarr, np.real(ptFJarr), 'r+', markersize=6)
        ax3.plot(karr, np.imag(Parr), lw=0.8, label='Imag') 
        ax3.plot(ptKarr, np.imag(ptFJarr), 'r+', markersize=6)
    else:
        if RorI:
            ax3.plot(karr, np.real(Parr), lw=0.8, label='Real') 
            ax3.plot(ptKarr, np.real(ptFJarr), 'r+', markersize=6)
        else:
            ax3.plot(karr, np.imag(Parr), lw=0.8, label='Imag') 
            ax3.plot(ptKarr, np.imag(ptFJarr), 'r+', markersize=6)
    

    ax3.set_title(rf'$\sum_k$ {FJname}')
    ax3.set_xlabel("k /$km^{-1}$")
    ax3.grid()
    ax3.legend(loc='lower left')

    return fig, (ax1, ax2, ax3)






def _dynamic_geometry_options(
    *, dist, azimuth, deprcv, src_fault, src_fault_size, rcv_points, rcv_fault, rcv_fault_size, nthreads,
) -> List[str]:
    """
    校验动态合成的接收几何，并格式化有限断层、任意点及线程参数

    :param    dist:            极坐标接收点的震中距
    :param    azimuth:         极坐标接收点的方位角
    :param    deprcv:          极坐标接收点的深度
    :param    src_fault:       震源断层文件
    :param    src_fault_size:  震源断层的剖分尺寸
    :param    rcv_points:      任意接收点文件
    :param    rcv_fault:       接收断层文件
    :param    rcv_fault_size:  接收断层的剖分尺寸
    :param    nthreads:        子源计算线程数
    :return:                  可直接传给 CLI 的选项列表
    """
    options = []
    explicit = sum((rcv_points is not None, rcv_fault is not None))
    if explicit > 1 or (explicit and (dist is not None or azimuth is not None)):
        raise ValueError("Receiver geometry options are mutually exclusive.")
    if (rcv_points is not None or rcv_fault is not None) and deprcv is not None:
        raise ValueError("rcv_points/rcv_fault supply their own depths; omit deprcv.")

    if not explicit and azimuth is None:
        raise ValueError("A polar receiver requires azimuth.")
    # 将有限断层或任意点参数直接加入命令
    for flag, path, size in (("C", src_fault, src_fault_size), ("U", rcv_fault, rcv_fault_size)):
        if size is not None and path is None:
            raise ValueError("Subdivision size requires the corresponding fault file.")
        if path is not None:
            option = f"-{flag}{Path(path)}"
            if size is not None:
                if len(size) != 2 or not (all(value == 0 for value in size) or all(value > 0 for value in size)):
                    raise ValueError("Fault size must contain two zeros or two positive values.")
                option += f"+i{format_float(size[0])}/{format_float(size[1])}"
            options.append(option)
    if rcv_points is not None:
        options.append(f"-Q{Path(rcv_points)}")

    if nthreads is not None:
        if isinstance(nthreads, bool) or int(nthreads) != nthreads or nthreads <= 0:
            raise ValueError("nthreads must be a positive integer.")
        options.append(f"-P{int(nthreads)}")
    return options


class _LambPhase(IntFlag):
    """与 C 的 GRT_LAMB_PHASE 保持一致的震相掩码"""
    P   = 1 << 0
    S   = 1 << 1
    R   = 1 << 2
    PP  = 1 << 3
    SS  = 1 << 4
    PS  = 1 << 5
    SP  = 1 << 6
    SPS = 1 << 7


def _prepare_lamb_phases(phases: Optional[Union[str, Sequence[str]]]) -> Optional[str]:
    """将 Python 侧的 Lamb 震相参数转换为逗号分隔的字符串"""
    if phases is None:
        return None
    if isinstance(phases, str):
        if phases == "":
            raise ValueError("phases should not be empty.")
        return phases
    try:
        phase_names = list(phases)
    except TypeError:
        raise TypeError("phases should be a string or a sequence of strings.") from None
    if not phase_names:
        raise ValueError("phases should contain at least one phase name.")
    if any(not isinstance(phase_name, str) for phase_name in phase_names):
        raise TypeError("each phase in phases should be a string.")
    phase_list = ",".join(phase_names)
    if phase_list == "":
        raise ValueError("phases should contain at least one non-empty phase name.")
    return phase_list


def _prepare_lamb_phase_mask(phases: Optional[Union[str, Sequence[str]]], supported: _LambPhase) -> int:
    """
    规范震相列表并在求解前解析为掩码

    :param    phases:     用户指定的震相列表，None 表示选择全部可用震相
    :param    supported:  当前求解器支持的震相掩码
    :return:             解析后的震相掩码
    """
    phase_list = _prepare_lamb_phases(phases)
    return C_grt_lamb_parse_phase_list(None if phase_list is None else phase_list.encode(), int(supported))


def lamb1(
    *, nu: float, tbar: np.ndarray, azimuth: float, cbar: Optional[float] = None,
    phases: Optional[Union[str, Sequence[str]]] = None,
):
    r"""
        solve the first-kind Lamb's problem using the generalized closed-form solution, see：

            张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社

        :param      nu:         Poisson ratio in (0, 0.5)
        :param      tbar:       dimensionless time :math:`\bar{t}=\dfrac{t}{T_S}=\dfrac{t}{r/\beta}=\dfrac{\beta t}{r}`,
                                where :math:`T_S=r/\beta` is the S-wave time scale and :math:`r` is the direct source-receiver distance
        :param      azimuth:    azimuth in degree, in ``[0, 360]``
        :param      cbar: dimensionless source velocity :math:`c/\beta` along positive ``x1``;
                                when specified, it must be positive and uses the Chapter 9
                                sub-Rayleigh moving-point-load solution for a vertical force source
        :param      phases: optional comma-separated phase selection or a non-empty sequence of phase names;
                                supported fixed-source phases are P, S and R

        :return:    Without ``cbar``, return the fixed-source Green function array ``G``
                    with shape ``(nt, 3, 3)``. When ``cbar`` is given, return the vertical-force
                    displacement array ``u`` with shape ``(nt, 3)``. Divide either result by
                    :math:`\pi^2 \mu r` to obtain the physical displacement, where :math:`\mu`
                    is the shear modulus and :math:`r` is the source-receiver distance.
    """

    nu = _prepare_lamb_scalar(nu, "nu")
    tbar = _prepare_lamb_time_series(tbar)
    azimuth = _prepare_lamb_scalar(azimuth, "azimuth")
    moving_source = cbar is not None
    cbar = 0.0 if cbar is None else _prepare_lamb_scalar(cbar, "cbar")
    if nu <= 0.0 or nu >= 0.5:
        raise ValueError("nu should be in (0, 0.5).")
    if azimuth < 0.0 or azimuth > 360.0:
        raise ValueError("azimuth should be in [0, 360].")
    if moving_source and cbar <= 0.0:
        raise ValueError("cbar should be positive when specified.")
    if cbar > 0.0 and abs(np.sin(np.deg2rad(azimuth))) <= 1e-8:
        raise ValueError("the moving-source closed-form solution requires azimuth off the x1 axis.")
    if moving_source and phases is not None:
        raise ValueError("phases is not supported with cbar.")

    phase_mask = _prepare_lamb_phase_mask(phases, _LambPhase.P | _LambPhase.S | _LambPhase.R)

    # 定义结果数组
    nt = len(tbar)
    u = np.zeros((nt, 3, 3), dtype=NPCT_REAL_TYPE)

    C_grt_solve_lamb1(
        nu,
        npct.as_ctypes(tbar),
        nt,
        azimuth,
        cbar,
        phase_mask,
        npct.as_ctypes(u.ravel()),
    )

    return u[:, :, 2].copy() if moving_source else u


def solve_lamb1(*args, **kwargs):
    """Legacy interface renamed to :func:`lamb1`; calling it raises an error."""
    raise RuntimeError("solve_lamb1() has been renamed to lamb1(); use lamb1() instead.")


def _prepare_lamb_scalar(value, name):
    try:
        value = np.asarray(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} should be a real scalar.") from None
    if value.ndim != 0:
        raise ValueError(f"{name} should be a real scalar.")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{name} should be a real scalar.") from None
    return value


def _prepare_lamb_time_series(tbar):
    try:
        tbar = np.asarray(tbar)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("tbar should be a one-dimensional sequence of real numbers.") from None
    if tbar.ndim != 1:
        raise ValueError("tbar should be a one-dimensional sequence of real numbers.")
    if np.iscomplexobj(tbar):
        raise ValueError("tbar should contain real values.")
    try:
        tbar = tbar.astype(NPCT_REAL_TYPE, copy=False)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("tbar should be a one-dimensional sequence of real numbers.") from None
    if tbar.size == 0:
        raise ValueError("tbar should not be empty.")
    if np.any(tbar < 0.0):
        raise ValueError("tbar should be nonnegative.")
    if tbar.size > 1 and np.any(np.diff(tbar) <= 0.0):
        raise ValueError("tbar should be strictly increasing.")
    return np.ascontiguousarray(tbar)


def _prepare_lamb_inputs(nu, tbar, R, azimuth):
    nu = _prepare_lamb_scalar(nu, "nu")
    tbar = _prepare_lamb_time_series(tbar)
    R = _prepare_lamb_scalar(R, "R")
    azimuth = _prepare_lamb_scalar(azimuth, "azimuth")
    if nu <= 0.0 or nu >= 0.5:
        raise ValueError("nu should be in (0, 0.5).")
    if R <= 0.0:
        raise ValueError("R should be positive.")
    if azimuth < 0.0 or azimuth > 360.0:
        raise ValueError("azimuth should be in [0, 360].")
    return nu, tbar, R, azimuth


def _prepare_lamb2_inputs(nu, tbar, R, depsrc, deprcv, azimuth):
    nu, tbar, R, azimuth = _prepare_lamb_inputs(nu, tbar, R, azimuth)
    src_given = depsrc is not None
    rcv_given = deprcv is not None
    if src_given == rcv_given:
        raise ValueError("exactly one of depsrc and deprcv should be set.")
    if src_given:
        depsrc = _prepare_lamb_scalar(depsrc, "depsrc")
        if depsrc <= 0.0:
            raise ValueError("depsrc should be strictly positive.")
        deprcv = 0.0
    else:
        deprcv = _prepare_lamb_scalar(deprcv, "deprcv")
        if deprcv <= 0.0:
            raise ValueError("deprcv should be strictly positive.")
        depsrc = 0.0
    return nu, tbar, R, depsrc, deprcv, azimuth


def lamb2(
    *, nu: float, tbar: np.ndarray, R: float,
    depsrc: Optional[float] = None, deprcv: Optional[float] = None, azimuth: float,
    phases: Optional[Union[str, Sequence[str]]] = None,
):
    r"""
        Solve the second-kind Lamb problem using the generalized closed-form solution.

        Exactly one of the source and receiver is on the free surface.
        Set ``depsrc`` to place the source underground and the receiver on the
        surface; set ``deprcv`` to place the receiver underground and the source
        on the surface, obtained from the buried-source solution by reciprocity.
        Exactly one of ``depsrc`` and ``deprcv`` may be set. ``R`` is the
        horizontal epicentral distance. The time array contains the dimensionless
        time :math:`\bar{t}=t/T_S=t/(r/\beta)=\beta t/r`, where
        :math:`T_S=r/\beta` and ``r = sqrt(R**2 + h**2)`` are the S-wave
        time scale and the straight source-receiver distance, and ``h`` is the
        nonzero depth. See:

            张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社

        :param      nu:           Poisson ratio in ``(0, 0.5)``; values within ``1e-3`` of either bound
                                    trigger a numerical warning
        :param      tbar:         dimensionless time :math:`\bar{t}=t/T_S=t/(r/\beta)=\beta t/r`,
                                    where :math:`T_S=r/\beta`
        :param      R:            positive horizontal epicentral distance from source to receiver;
                                    when ``R / sqrt(R**2 + h**2) <= 1e-3``,
                                    a numerical-stability warning is issued
        :param      depsrc:       strictly positive source depth with the receiver on the surface;
                                    mutually exclusive with ``deprcv``. Values below
                                    ``1e-3 * r`` trigger a numerical warning
        :param      deprcv:       strictly positive receiver depth with the source on the surface;
                                    mutually exclusive with ``depsrc``. Values below
                                    ``1e-3 * r`` trigger a numerical warning
        :param      azimuth:      azimuth in degree, from source to receiver, in ``[0, 360]``
        :param      phases:       optional comma-separated phase selection or a non-empty sequence of phase names;
                                    supported phases are P, S, SP and PS
        :return:    Four normalized arrays ``G, Gs, Gr, Grs``. ``G`` has shape
                    ``(nt, 3, 3)`` and the first-derivative arrays have shape
                    ``(nt, 3, 3, 3)``. ``Grs`` has shape ``(nt, 3, 3, 3, 3)``.
                    The index order is ``[time, coordinate, receiver_component,
                    source_component]`` for first derivatives and
                    ``[time, receiver_coordinate, source_coordinate, receiver_component,
                    source_component]`` for ``Grs``. To obtain the physical
                    solutions, divide ``G`` by :math:`\pi^2 \mu r`, divide ``Gs``
                    and ``Gr`` by :math:`\pi^2 \mu r^2`, and divide ``Grs`` by
                    :math:`\pi^2 \mu r^3`, where :math:`\mu` is the shear modulus
                    and :math:`r` is the source-receiver distance.
    """

    nu, tbar, R, depsrc, deprcv, azimuth = _prepare_lamb2_inputs(nu, tbar, R, depsrc, deprcv, azimuth)
    supported = _LambPhase.P | _LambPhase.S | (_LambPhase.SP if depsrc > 0 else _LambPhase.PS)
    phase_mask = _prepare_lamb_phase_mask(phases, supported)
    nt = len(tbar)
    G = np.zeros((nt, 3, 3), dtype=NPCT_REAL_TYPE)
    Gr = np.zeros((nt, 3, 3, 3), dtype=NPCT_REAL_TYPE)
    Gs = np.zeros((nt, 3, 3, 3), dtype=NPCT_REAL_TYPE)
    Grs = np.zeros((nt, 3, 3, 3, 3), dtype=NPCT_REAL_TYPE)

    C_grt_solve_lamb2(
        nu,
        npct.as_ctypes(tbar),
        nt,
        R,
        depsrc,
        deprcv,
        azimuth,
        phase_mask,
        npct.as_ctypes(G.ravel()),
        npct.as_ctypes(Gs.ravel()),
        npct.as_ctypes(Gr.ravel()),
        npct.as_ctypes(Grs.ravel()),
    )

    return G, Gs, Gr, Grs


def _prepare_lamb3_inputs(nu, tbar, R, depsrc, deprcv, azimuth):
    nu, tbar, R, azimuth = _prepare_lamb_inputs(nu, tbar, R, azimuth)
    depsrc = _prepare_lamb_scalar(depsrc, "depsrc")
    deprcv = _prepare_lamb_scalar(deprcv, "deprcv")
    if depsrc <= 0.0 or deprcv <= 0.0:
        raise ValueError("depsrc and deprcv should be strictly positive.")
    return nu, tbar, R, depsrc, deprcv, azimuth


def lamb3(
    *, nu: float, tbar: np.ndarray, R: float, depsrc: float, deprcv: float, azimuth: float,
    phases: Optional[Union[str, Sequence[str]]] = None,
):
    r"""
        Solve the third-kind Lamb problem using the generalized closed-form solution.

        Both the source and receiver are inside the halfspace. ``R`` is their
        horizontal epicentral distance, and the two depths use the same length unit.
        The time array contains the dimensionless time
        :math:`\bar{t}=t/T_S=t/(r/\beta)=\beta t/r`, where
        :math:`T_S=r/\beta` and ``r = sqrt(R**2 + (depsrc - deprcv)**2)``
        are the S-wave time scale and the straight source-receiver distance,
        respectively.

        See:

            张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社

        :param      nu:           Poisson ratio in ``(0, 0.5)``; values within ``1e-3`` of either bound
                                    trigger a numerical warning
        :param      tbar:         dimensionless time :math:`\bar{t}=t/T_S=t/(r/\beta)=\beta t/r`,
                                    where :math:`T_S=r/\beta`
        :param      R:            positive horizontal source-receiver distance; when
                                    ``R / sqrt(R**2 + (depsrc + deprcv)**2) <= 1e-2``,
                                    a numerical-stability warning is issued
        :param      depsrc:       strictly positive source depth; values below ``1e-3 * r``
                                    trigger a numerical warning
        :param      deprcv:       strictly positive receiver depth; values below ``1e-3 * r``
                                    trigger a numerical warning
        :param      azimuth:      azimuth in degree, from source to receiver, in ``[0, 360]``
        :param      phases:       optional comma-separated phase selection or a non-empty sequence of phase names;
                                    supported phases are P, S, PP, SS, PS, SP and sPs
        :return:    Four normalized arrays ``G, Gs, Gr, Grs``. ``G`` has shape
                    ``(nt, 3, 3)`` and the first-derivative arrays have shape
                    ``(nt, 3, 3, 3)``. ``Grs`` has shape ``(nt, 3, 3, 3, 3)``.
                    The index order is ``[time, coordinate, receiver_component,
                    source_component]`` for first derivatives and
                    ``[time, receiver_coordinate, source_coordinate, receiver_component,
                    source_component]`` for ``Grs``. To obtain the physical
                    solutions, divide ``G`` by :math:`\pi^2 \mu r`, divide ``Gs``
                    and ``Gr`` by :math:`\pi^2 \mu r^2`, and divide ``Grs`` by
                    :math:`\pi^2 \mu r^3`, where :math:`\mu` is the shear modulus
                    and :math:`r` is the source-receiver distance.
    """

    nu, tbar, R, depsrc, deprcv, azimuth = _prepare_lamb3_inputs(nu, tbar, R, depsrc, deprcv, azimuth)
    supported = _LambPhase.P | _LambPhase.S | _LambPhase.PP | _LambPhase.SS | _LambPhase.PS | _LambPhase.SP | _LambPhase.SPS
    phase_mask = _prepare_lamb_phase_mask(phases, supported)
    nt = len(tbar)
    G = np.zeros((nt, 3, 3), dtype=NPCT_REAL_TYPE)
    Gr = np.zeros((nt, 3, 3, 3), dtype=NPCT_REAL_TYPE)
    Gs = np.zeros((nt, 3, 3, 3), dtype=NPCT_REAL_TYPE)
    Grs = np.zeros((nt, 3, 3, 3, 3), dtype=NPCT_REAL_TYPE)
    C_grt_solve_lamb3(
        nu,
        npct.as_ctypes(tbar),
        nt,
        R,
        depsrc,
        deprcv,
        azimuth,
        phase_mask,
        npct.as_ctypes(G.ravel()),
        npct.as_ctypes(Gs.ravel()),
        npct.as_ctypes(Gr.ravel()),
        npct.as_ctypes(Grs.ravel()),
    )
    return G, Gs, Gr, Grs


def lamb(
    *,
    modelparams: Sequence[float],
    depsrc: Optional[float] = None,
    deprcv: Optional[float] = None,
    dist: Optional[float] = None,
    nt: int,
    dt: float,
    azimuth: Optional[float] = None,
    output_path: PathLike,
    scale: Optional[float] = None,
    scale_with_mu: bool = False,
    strike: Optional[float] = None,
    dip: Optional[float] = None,
    rake: Optional[float] = None,
    force: Optional[Sequence[float]] = None,
    moment_tensor: Optional[Sequence[float]] = None,
    src_fault: Optional[PathLike] = None,
    src_fault_size: Optional[Sequence[float]] = None,
    rcv_fault: Optional[PathLike] = None,
    rcv_fault_size: Optional[Sequence[float]] = None,
    rcv_points: Optional[PathLike] = None,
    nthreads: Optional[int] = None,
    time_function: Optional[str] = None,
    integrate_order: Optional[int] = None,
    differentiate_order: Optional[int] = None,
    phases: Optional[Union[str, Sequence[str]]] = None,
    delayT0: float = 0.0,
    delayV0: float = 0.0,
    ref_first_p: bool = False,
    zne: Optional[bool] = None,
    calc_upar: bool = False,
    print_log: bool = True,
) -> None:
    r"""
    Synthesize dynamic displacement with the physical Lamb closed-form solution.

    A call defines five groups of information:

    * Medium and sampling: ``modelparams`` gives ``(vp, vs, rho)`` for the homogeneous
      half-space; ``nt`` and ``dt`` set the output sampling. No Green-function library is needed.
    * Source location: a point source is at the horizontal origin, with depth ``depsrc``.
      ``src_fault`` supplies the locations of finite sources.
    * Receiver locations: choose polar coordinates (``dist``, ``azimuth``, ``deprcv``),
      a point file (``rcv_points``), or receiver faults (``rcv_fault``).
      Use one receiver mode per call; files supply their own depths.
    * Source mechanism and strength: a point source requires ``scale``.
      Choose ``strike``/``dip``/``rake`` for a double-couple, omit ``rake`` for a tensile crack,
      or use ``force`` or ``moment_tensor`` instead. Leaving all mechanism parameters unset
      selects an explosion. ``src_fault`` supplies both location and mechanism/strength;
      omit ``depsrc``, ``scale`` and all point-source mechanism parameters in this mode.
    * Output: ``output_path`` sets the SAC directory; ``zne`` selects ZNE
      and ``calc_upar`` adds spatial derivatives. Finite sources always output ZNE.

    Point sources require ``depsrc``; polar receivers require all three of ``dist``,
    ``azimuth`` and ``deprcv``. All coordinates share one horizontal origin;
    for finite sources, ``dist``/``azimuth`` locate the receiver relative to that origin.
    Rectangular sources require ``src_fault_size``; set ``(0, 0)`` to use one center
    point source weighted by the whole fault area. Receiver subdivision is optional;
    omit ``rcv_fault_size`` or set ``(0, 0)`` to use one center point per fault.

    Optional ``time_function``, ``integrate_order`` and ``differentiate_order`` control
    the time dependence. ``delayT0``, ``delayV0`` and ``ref_first_p`` set the output start;
    ``phases`` selects phases and ``nthreads`` sets the source thread count.
    The output sample count is the maximum of ``nt`` and the longest source time
    function including rupture delay, shared by all receivers. Convolution is linear.
    The source and receiver depths select the appropriate Lamb solution.
    When both depths are zero, only ``force`` is supported. ``calc_upar`` is ignored
    for a surface force with all receivers on the surface; a mixture of surface
    and buried receivers with ``calc_upar=True`` raises an error.
    All arguments are keyword-only. Results are displacements in cm with Z upward,
    R radial outward and T clockwise from R by default. Polar receivers write directly
    under the output path; receiver files always use indexed subdirectories, even for
    one point. Finite-source ``sig.sac`` stores the total scalar moment rate
    only when a row or global time function is explicitly specified.
    A global ``time_function`` ignores all row-end contents in the source fault file.
    Results are written to files; read SAC waveforms explicitly with :func:`obspy.read`.

    Each receiver records its earliest P and S arrivals over nonzero source points,
    including sampled explicit rupture delays. With multiple source points, later
    phase picks are left undefined. Arrivals and output starts are determined after
    fault subdivision. The internal start is fixed near the origin time; the output window
    sets only the solve end time. Convolution and time integration/differentiation retain
    the full response history, then the result is shifted by rupture delay and cropped.

    :param    modelparams:      Homogeneous half-space parameters ``(vp, vs, rho)``;
                                velocities are in km/s and density is in g/cm^3
    :param    depsrc:           Point-source depth in km. Required for a point source; omit for
                                source faults.
    :param    deprcv:           Polar receiver depth in km. Required for polar receivers; omit
                                for receiver files.
    :param    dist:             Positive polar receiver distance from the horizontal origin in
                                km. Required for polar receivers; omit for receiver files.
    :param    nt:               Minimum number of time samples; extended to cover the longest
                                source time function including rupture delay, equally for all receivers
    :param    dt:               Time-sample interval in s
    :param    azimuth:          Polar receiver azimuth in degrees clockwise from north. Required
                                for polar receivers; omit for receiver files.
    :param    output_path:      Output directory for SAC files
    :param    scale:            Source scaling factor
    :param    scale_with_mu:    Whether to multiply ``scale`` by the source-layer
                                shear modulus
    :param    strike:           Fault strike in degrees
    :param    dip:              Fault dip in degrees
    :param    rake:             Slip rake in degrees
    :param    force:            Single-force coefficients ``(fN, fE, fZ)``
    :param    moment_tensor:    Moment-tensor coefficients
                                ``(Mxx, Mxy, Mxz, Myy, Myz, Mzz)``
    :param    src_fault:        Coulomb source file supplying source positions, mechanisms and
                                slip/potency. Append a complete ``-D`` option to each row
                                to specify its rupture process.
    :param    src_fault_size:   Along-strike/dip subdivision sizes (dL, dW) in km,
                                required for rectangular sources. Both values must be positive
                                or zero; ``(0, 0)`` disables further subdivision.
                                Point Kode records remain single points.
    :param    rcv_points:       ASCII receiver file: north east depth in km, optionally followed
                                by strike dip rake in degrees.
    :param    rcv_fault:        Coulomb receiver file; only Kode=100 is supported and
                                slip magnitude is ignored.
                                An exact ``rake`` header preserves the angle even at zero slip;
                                otherwise the slip columns define direction.
    :param    rcv_fault_size:   Along-strike/dip receiver subdivision sizes (dL, dW) in km.
                                Both values must be positive or zero.
                                Omit or set ``(0, 0)`` to use one center point per fault.
    :param    nthreads:         Positive OpenMP source-subfault thread count.
    :param    time_function:    Time-function parameters passed to ``grt``, without the ``-D`` prefix.
                                Supported forms are ``i`` (impulse), ``p/t0``, ``t/t1/t2/t3``,
                                ``c/t1/t2`` (asymmetric cosine) or ``0/file``, with area normalization.
                                Custom files support one amplitude column or two columns: time (s) and
                                amplitude. Two-column data are linearly interpolated to ``dt``.
                                Custom time functions are automatically area-normalized.
                                For ``t/t1/t2/t3``, the parameters are nonnegative rise, plateau
                                and fall durations in seconds, with a positive total duration.
                                For ``c/t1/t2``, ``t1`` and ``t2`` are positive rise and fall durations in seconds.
                                ``r/f0`` additionally accepts a signed Ricker convolution wavelet,
                                with ``f0`` in Hz and analytic peak amplitude 1, without area normalization.
                                It is not a unit-slip source process.
                                Append ``+d<delay>`` for a delay in seconds, e.g. ``p/1.3+d0.4``.
    :param    integrate_order:  Number of time integrations
    :param    differentiate_order: Number of time differentiations
    :param    phases:            Optional comma-separated phase selection or a non-empty sequence of phase names;
                                the available names depend on the selected Lamb problem
    :param    delayT0:          Time delay at zero distance in s
    :param    delayV0:          Reference velocity in km/s. The reference distance is the shortest
                                straight distance from a nonzero source point to the receiver.
    :param    ref_first_p:      Whether to reference the output start to the earliest P arrival,
                                including sampled rupture delays. Requires negative delayT0.
    :param    zne:              If true, output ZNE instead of ZRT. Finite sources always use
                                ZNE.
    :param    calc_upar:        Whether to output spatial displacement derivatives
    :param    print_log:        Whether to print regular ``grt`` output
    """
    vp, vs, rho = modelparams
    output = Path(output_path)
    command = [
        "lamb",
        f"-H{format_float(vp)}/{format_float(vs)}/{format_float(rho)}",
        f"-N{nt}/{format_float(dt)}",
        f"-O{output}",
    ]
    command.extend(_dynamic_geometry_options(
        dist=dist, azimuth=azimuth, deprcv=deprcv, src_fault=src_fault, src_fault_size=src_fault_size,
        rcv_points=rcv_points, rcv_fault=rcv_fault, rcv_fault_size=rcv_fault_size, nthreads=nthreads,
    ))
    if src_fault is not None:
        if any(value is not None for value in (depsrc, scale, strike, dip, rake, force, moment_tensor)) or scale_with_mu:
            raise ValueError("src_fault is mutually exclusive with point-source parameters.")
        if zne is False:
            warnings.warn("Finite sources always output ZNE; zne=False is ignored.", stacklevel=2)
        zne = True
    elif scale is None:
        raise ValueError("Point sources require scale.")
    for flag, value in (("R", dist), ("Ds", depsrc), ("Dr", deprcv), ("A", azimuth)):
        if value is not None:
            command.append(f"-{flag}{format_float(value)}")
    if scale is not None:
        command.append(f"-S{'u' if scale_with_mu else ''}{format_float(scale)}")
    phase_list = _prepare_lamb_phases(phases)
    if phase_list is not None:
        command.append(f"-L{phase_list}")

    # 点源机制由单力、矩张量或断层角度三种互斥输入确定
    has_mechanism = strike is not None or dip is not None or rake is not None
    if force is not None:
        if has_mechanism or moment_tensor is not None:
            raise ValueError("force is mutually exclusive with strike/dip/rake and moment_tensor.")
        command.append("-F" + "/".join(format_float(value) for value in force))
    elif moment_tensor is not None:
        if has_mechanism:
            raise ValueError("moment_tensor is mutually exclusive with strike/dip/rake and force.")
        command.append("-T" + "/".join(format_float(value) for value in moment_tensor))
    elif has_mechanism:
        if strike is None or dip is None:
            raise ValueError("strike and dip must be supplied together.")
        source = f"-M{format_float(strike)}/{format_float(dip)}"
        if rake is not None:
            source += f"/{format_float(rake)}"
        command.append(source)

    if time_function is not None:
        command.append(f"-D{time_function}")
    if ref_first_p:
        command.append(f"-Ep{format_float(delayT0)}")
    else:
        command.append(f"-E{format_float(delayT0)}/{format_float(delayV0)}")
    if integrate_order is not None:
        command.append(f"-I{integrate_order}")
    if differentiate_order is not None:
        command.append(f"-J{differentiate_order}")
    if zne:
        command.append("-n")
    if calc_upar:
        command.append("-e")

    output.mkdir(parents=True, exist_ok=True)
    run_grt(command, print_log=print_log)
    return None


# ======================================================================================================
#                                           面波辅助模块
# ======================================================================================================


def _is_rayl(data: dict) -> bool:
    """从 attributes 或本征函数列数判断是否为 Rayleigh 波"""
    val = data.get("attributes", {}).get("isRayl")
    if val is not None:
        return bool(np.asarray(val).reshape(-1)[0])
    eigfn = data.get("eigfn")
    if eigfn is None:
        return True
    if isinstance(eigfn, dict):
        if not eigfn:
            return True
        eigfn = next(iter(eigfn.values()))
    return int(np.asarray(eigfn).shape[-1]) >= 8


def _as_modal_panels(data: dict, title: Optional[str] = None, required: Sequence[str] = ("freq",)) -> list[tuple[str, dict]]:
    """将单个结果或 {标题: 结果} 统一成 [(title, data), ...]"""
    def wave_title(data: dict) -> str:
        attrs = data.get("attributes", {})
        if "isRayl" in attrs:
            return "Rayleigh" if _is_rayl(data) else "Love"
        return ""

    if not isinstance(data, dict) or not data:
        raise TypeError("data must be a non-empty dict.")
    if all(key in data for key in required):
        if not title:
            title = wave_title(data)
        return [(title, data)]

    panels = []
    for key, val in data.items():
        if not isinstance(val, dict) or not all(name in val for name in required):
            raise ValueError("Mapping values must be dictionaries returned by the corresponding read function.")
        panels.append((str(key), val))
    return panels


def _prepare_modal_axes(panels: list[tuple[str, dict]], ax: Optional[Axes], figsize: tuple[float, float]):
    """创建或复用坐标轴，返回 (fig, axes_tuple)"""
    n = len(panels)
    if ax is not None:
        if n != 1:
            raise ValueError("ax can only be used with a single dataset.")
        return ax.figure, (ax,)

    fig, axs = plt.subplots(1, n, figsize=figsize, layout="constrained")
    if n == 1:
        axs = (axs,)
    else:
        axs = tuple(axs)
    return fig, axs


def _maybe_savefig(fig: Figure, outpath: Optional[PathLike]):
    """若给定路径则保存图片"""
    if outpath is not None:
        fig.savefig(outpath)


def _read_modal(path: PathLike, kind: str, required: Sequence[str], any_of: Sequence[str] = (), forbidden: Sequence[str] = ()) -> dict:
    """读取面波 NetCDF，并将沿 freqmode 展平的变量按阶号分组"""
    def group_by_mode(raw: dict, attributes: dict) -> dict:
        freq = np.asarray(raw["freq"])
        cnum = np.asarray(raw["cnum"])
        mode = np.asarray(raw["mode"])
        if freq.size != cnum.size:
            raise ValueError("freq and cnum must have the same length.")

        # 第 iw 个频率对应 mode 的前 cnum[iw] 个阶，频散点沿 freqmode 展平
        freq_of: dict[int, list] = {}
        idx_of: dict[int, list] = {}
        k = 0
        for iw, n in enumerate(cnum):
            for ic in range(int(n)):
                im = int(mode[ic])
                freq_of.setdefault(im, []).append(float(freq[iw]))
                idx_of.setdefault(im, []).append(k)
                k += 1
        nfm = k

        result = {
            "mode": np.array(sorted(idx_of), dtype=int),
            "freq": {im: np.asarray(freq_of[im], dtype=float) for im in idx_of},
        }
        freqmode_vars = ("c", "u", "ciref", "eigfn", "csens", "usens", "egyint")
        skip = {"freq", "cnum", "mode"}
        for name, arr in raw.items():
            if name in skip:
                continue
            arr = np.asarray(arr)
            if name in freqmode_vars:
                if arr.ndim < 1 or arr.shape[0] != nfm:
                    raise ValueError(f"{name} length does not match sum(cnum).")
                result[name] = {im: arr[np.asarray(idx_of[im])] for im in idx_of}
            else:
                result[name] = arr
        result["attributes"] = dict(attributes)
        return result

    path = str(path)
    nc = read_nc(path)
    raw = {name: info["data"] for name, info in nc["variables"].items()}
    missing = [name for name in required if name not in raw]
    if missing:
        raise ValueError(f"{path} is not a {kind} NetCDF file (missing {', '.join(missing)}).")
    if any_of and not any(name in raw for name in any_of):
        raise ValueError(f"{path} is not a {kind} NetCDF file (missing {'/'.join(any_of)}).")
    hit = [name for name in forbidden if name in raw]
    if hit:
        raise ValueError(f"{path} is not a {kind} NetCDF file.")
    return group_by_mode(raw, nc["attributes"])


def read_dispersion(path: PathLike) -> dict:
    """
    Read a phase- or group-velocity dispersion NetCDF file.

    Dispersion points are grouped by mode number. Typical keys:

    * ``freq`` - mapping from mode number to frequency (Hz)
    * ``c`` - mapping from mode number to phase velocity (km/s)
    * ``u`` - mapping from mode number to group velocity when present
    * ``mode`` - array of mode numbers present in the file
    * ``attributes`` - global NetCDF attributes

    :param    path:               Path to the dispersion NetCDF file.

    :return: A dictionary of dispersion results grouped by mode.
    """
    return _read_modal(
        path, "dispersion",
        required=("freq", "cnum", "mode"),
        any_of=("c", "u"),
        forbidden=("eigfn", "csens", "usens", "egyint"),
    )


def read_eigenfunction(path: PathLike) -> dict:
    """
    Read an eigenfunction NetCDF file produced by ``grt eigenfn -W``.

    Arrays along the flattened ``freqmode`` axis are grouped by mode.
    Typical keys:

    * ``freq`` - mapping from mode number to frequency (Hz)
    * ``c`` - mapping from mode number to phase velocity (km/s)
    * ``eigfn`` - mapping from mode number to array of shape ``(nf, nz, nw)``
    * ``z`` - depth array (km)
    * ``mode`` - array of mode numbers present in the file
    * ``attributes`` - global NetCDF attributes

    :param    path:               Path to the eigenfunction NetCDF file.

    :return: A dictionary of eigenfunction results grouped by mode.
    """
    return _read_modal(
        path, "eigenfunction",
        required=("freq", "cnum", "mode", "z", "eigfn"),
    )


def read_energy_integral(path: PathLike) -> dict:
    """
    Read an energy-integral NetCDF file produced by ``grt eigenfn -K+x``.

    Arrays along the flattened ``freqmode`` axis are grouped by mode.
    Typical keys:

    * ``freq`` - mapping from mode number to frequency (Hz)
    * ``c`` - mapping from mode number to phase velocity (km/s)
    * ``egyint`` - mapping from mode number to array of shape ``(nf, 10)``
    * ``mode`` - array of mode numbers present in the file
    * ``attributes`` - global NetCDF attributes

    :param    path:               Path to the energy-integral NetCDF file.

    :return: A dictionary of energy-integral results grouped by mode.
    """
    return _read_modal(
        path, "energy-integral",
        required=("freq", "cnum", "mode", "egyint"),
        forbidden=("eigfn", "csens", "usens"),
    )


def read_sensitivity(path: PathLike) -> dict:
    """
    Read a phase- or group-velocity sensitivity-kernel NetCDF file.

    Arrays along the flattened ``freqmode`` axis are grouped by mode.
    Typical keys:

    * ``freq`` - mapping from mode number to frequency (Hz)
    * ``csens`` / ``usens`` - mapping from mode number to array of shape ``(nf, nz, 3)``
    * ``z`` - depth array (km)
    * ``mode`` - array of mode numbers present in the file
    * ``attributes`` - global NetCDF attributes

    :param    path:               Path to the sensitivity NetCDF file.

    :return: A dictionary of sensitivity results grouped by mode.
    """
    return _read_modal(
        path, "sensitivity",
        required=("freq", "cnum", "mode", "z"),
        any_of=("csens", "usens"),
    )


def read_secfunc(path: PathLike) -> list[dict]:
    """
    Read a secular-function text file produced by ``grt eigenv -X``.

    Each comment line starting with ``#`` begins a new secular function
    (typically one layer ``iref``). Data lines have three columns: phase
    velocity (km/s), real part and imaginary part.

    :param    path:               Path to the secular-function text file.

    :return: A list of dictionaries. Each item contains ``c``, ``real``,
             ``imag``, and header fields such as ``iref`` and ``freq``.
    """
    def parse_header(line: str) -> dict:
        attrs = {}
        for part in line.lstrip("#").split(","):
            if "=" not in part:
                continue
            key, val = part.split("=", 1)
            key = key.strip()
            val = val.strip()
            if key == "iref":
                attrs["iref"] = int(float(val))
            elif key == "f":
                attrs["freq"] = float(val)
            elif key in ("tol", "dc", "c1", "c2"):
                attrs[key] = float(val)
        return attrs

    path = str(path)
    if not Path(path).is_file():
        raise FileNotFoundError(f"Secular-function file does not exist: {path}")

    records = []
    header = {}
    c_list, re_list, im_list = [], [], []

    def flush():
        if not c_list:
            return
        rec = {
            "c":    np.asarray(c_list, dtype=float),
            "real": np.asarray(re_list, dtype=float),
            "imag": np.asarray(im_list, dtype=float),
        }
        rec.update(header)
        records.append(rec)
        c_list.clear()
        re_list.clear()
        im_list.clear()

    with open(path, "r") as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            # 以 # 开头的注释行分隔不同层的久期函数
            if line.startswith("#"):
                flush()
                header = parse_header(line)
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(f"{path} is not a secular-function file (expect 3 columns).")
            c_list.append(float(parts[0]))
            re_list.append(float(parts[1]))
            im_list.append(float(parts[2]))
    flush()

    if not records:
        raise ValueError(f"{path} contains no secular-function data.")
    return records


def plot_dispersion(
    data: dict,
    *,
    varname: Optional[str] = None,
    title: Optional[str] = None,
    ax: Optional[Axes] = None,
    outpath: Optional[PathLike] = None,
):
    """
    Plot phase- or group-velocity dispersion curves.

    :param    data:               Dictionary from :func:`read_dispersion <pygrt.utils.read_dispersion>`,
                                  or a mapping from panel title to such dictionaries
                                  (drawn side by side).
    :param    varname:            Variable to plot, ``c`` for phase velocity or ``u`` for group
                                  velocity. Default is ``u`` when present, otherwise ``c``.
    :param    title:              Panel title when ``data`` is a single dictionary.
    :param    ax:                 Existing matplotlib Axes. Only valid for a single dataset.
    :param    outpath:            If given, save the figure to this path.

    :return:
            - **fig** -                        matplotlib.Figure object
            - **axs** -                        tuple of matplotlib.Axes
    """
    def resolve_varname(data: dict, varname: Optional[str]) -> str:
        if varname is None:
            varname = "u" if "u" in data else "c"
        if varname not in ("c", "u"):
            raise ValueError("varname must be 'c' or 'u'.")
        if varname not in data:
            raise KeyError(f"variable {varname!r} is not in the dispersion dictionary.")
        return varname

    def plot_on_ax(data: dict, ax: Axes, varname: str, title: str):
        for im in sorted(data[varname]):
            ax.plot(data["freq"][im], data[varname][im], lw=0.6, marker="o", ms=0, c="b")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Group velocity (km/s)" if varname == "u" else "Phase velocity (km/s)")
        if title:
            ax.set_title(title)
        ax.set_xmargin(0)
        ax.set_xlim(xmin=0)

    panels = _as_modal_panels(data, title)
    varname = resolve_varname(panels[0][1], varname)
    fig, axs = _prepare_modal_axes(panels, ax, figsize=(5.0 * max(len(panels), 1), 4.0))
    for (panel_title, panel_data), axi in zip(panels, axs):
        plot_on_ax(panel_data, axi, varname, panel_title)
    _maybe_savefig(fig, outpath)
    return fig, axs


def plot_eigenfunction(
    data: dict,
    *,
    scale: float = 2.0,
    ifreq: int = 0,
    title: Optional[str] = None,
    ax: Optional[Axes] = None,
    outpath: Optional[PathLike] = None,
):
    """
    Plot eigenfunctions of different modes at one frequency.

    Rayleigh-wave panels show the vertical component as a solid line and the
    radial component as a dashed line. Love-wave panels show the transverse
    component. Each mode is offset horizontally by its mode number.

    :param    data:               Dictionary from
                                  :func:`read_eigenfunction <pygrt.utils.read_eigenfunction>`,
                                  or a mapping from panel title to such dictionaries
                                  (drawn side by side).
    :param    scale:              Amplitude scale of the plotted eigenfunctions.
    :param    ifreq:              Index among frequencies that have at least one mode.
                                  Default is the first frequency.
    :param    title:              Panel title when ``data`` is a single dictionary.
    :param    ax:                 Existing matplotlib Axes. Only valid for a single dataset.
    :param    outpath:            If given, save the figure to this path.

    :return:
            - **fig** -                        matplotlib.Figure object
            - **axs** -                        tuple of matplotlib.Axes
    """
    def plot_on_ax(data: dict, ax: Axes, title: str):
        zarr = np.asarray(data["z"])
        eigfn = data["eigfn"]
        freq = data["freq"]
        all_f = np.unique(np.concatenate([np.asarray(v, dtype=float) for v in freq.values()]))
        ifreq_i = ifreq + all_f.size if ifreq < 0 else ifreq
        if ifreq_i < 0 or ifreq_i >= all_f.size:
            raise IndexError(f"ifreq={ifreq} is out of range for nf={all_f.size}.")
        target = float(all_f[ifreq_i])

        is_rayl = _is_rayl(data)
        idx_vert, idx_horz = 2, 0
        solidkwargs = dict(c="k", ls="-", lw=0.5)
        dashkwargs = dict(c="k", ls="--", lw=0.5)

        for im in sorted(eigfn):
            loc = np.flatnonzero(np.isclose(np.asarray(freq[im], dtype=float), target))
            if loc.size == 0:
                continue
            fn = np.asarray(eigfn[im][loc[0]])
            if is_rayl:
                trace = np.array(fn[:, idx_vert], dtype=float, copy=True)
                norm = np.max(np.abs(trace))
                if norm == 0.0:
                    norm = 1.0
                trace /= norm * scale
                ax.plot(im + trace, zarr, **solidkwargs)
                trace = np.array(fn[:, idx_horz], dtype=float, copy=True)
                trace /= norm
                ax.plot(im + trace, zarr, **dashkwargs)
            else:
                trace = np.array(fn[:, idx_horz], dtype=float, copy=True)
                norm = np.max(np.abs(trace))
                if norm == 0.0:
                    norm = 1.0
                trace /= norm * scale
                ax.plot(im + trace, zarr, **solidkwargs)

        ax.yaxis.set_inverted(True)
        ax.grid()
        ax.set_ymargin(0)
        ax.set_ylabel("Depth (km)")
        ax.set_xlabel("Order of modes")
        if title:
            ax.set_title(title)

        if is_rayl:
            ax.legend(
                [Line2D([], [], **solidkwargs), Line2D([], [], **dashkwargs)],
                ["Vertical", "Radial"],
                loc="lower right", ncol=1,
            )
        else:
            ax.legend([Line2D([], [], **solidkwargs)], ["Transverse"], loc="lower right", ncol=1)

    if scale <= 0.0:
        raise ValueError("scale must be positive.")
    panels = _as_modal_panels(data, title, required=("freq", "eigfn"))
    fig, axs = _prepare_modal_axes(panels, ax, figsize=(5.0 * max(len(panels), 1), 4.0))
    for (panel_title, panel_data), axi in zip(panels, axs):
        plot_on_ax(panel_data, axi, panel_title)
    _maybe_savefig(fig, outpath)
    return fig, axs


def plot_sensitivity(
    data: dict,
    *,
    modes: Union[int, Sequence[int], None] = 0,
    title: Optional[str] = None,
    outpath: Optional[PathLike] = None,
):
    r"""
    Plot dimensionless phase- or group-velocity sensitivity kernels.

    Three panels correspond to :math:`\alpha`, :math:`\beta` and :math:`\rho`.
    By default only the fundamental mode (``modes=0``) is drawn; pass ``None``
    to draw all modes.

    :param    data:               Dictionary from
                                  :func:`read_sensitivity <pygrt.utils.read_sensitivity>`.
    :param    modes:              Mode number, a sequence of mode numbers, or ``None``
                                  for all modes. Default is the fundamental mode.
    :param    title:              Figure suptitle.
    :param    outpath:            If given, save the figure to this path.

    :return:
            - **fig** -                        matplotlib.Figure object
            - **axs** -                        tuple of matplotlib.Axes
    """
    def normalize_modes(modes) -> Optional[set[int]]:
        if modes is None:
            return None
        if isinstance(modes, (int, np.integer)):
            return {int(modes)}
        return {int(m) for m in modes}

    def plot_component(idx: int, ax: Axes):
        zs = np.asarray(data["z"])
        hLst = []
        for im in sorted(sens):
            if mode_set is not None and im not in mode_set:
                continue
            arr = np.asarray(sens[im])
            farr = np.asarray(data["freq"][im])
            for i in range(arr.shape[0]):
                freq = float(farr[i])
                if multi_mode:
                    label = rf"$n={im},\ f={freq:.1f}\ \mathrm{{Hz}}$"
                else:
                    label = f"$f={freq:.1f} Hz$"
                h, = ax.plot(arr[i, :, idx], zs, lw=0.3, marker="o", ms=2, label=label)
                hLst.append(h)
        ax.set_ymargin(0)
        ax.yaxis.set_inverted(True)
        return hLst

    if not isinstance(data, dict) or "freq" not in data:
        raise TypeError("data must be a dictionary returned by read_sensitivity().")
    if "csens" in data:
        char, sname = "c", "csens"
    elif "usens" in data:
        char, sname = "U", "usens"
    else:
        raise KeyError("dictionary has neither csens nor usens.")

    sens = data[sname]
    mode_set = normalize_modes(modes)
    multi_mode = mode_set is None or len(mode_set) != 1
    symbols = [r"\alpha", r"\beta", r"\rho"]

    with plt.rc_context({"mathtext.fontset": "cm"}):
        fig, axs = plt.subplots(1, 3, figsize=(5, 6), layout="constrained", sharey=True)
        hLst = []
        for i, axi in enumerate(axs):
            handles = plot_component(i, axi)
            axi.set_xlabel(rf"$\dfrac{{{symbols[i]}}}{{{char}}} \dfrac{{\partial {char}}}{{\partial {symbols[i]}}}$", fontsize=14)
            if i == 1:
                hLst = handles
            if i == 0:
                axi.set_ylabel("Depth (km)")
        if hLst:
            axs[1].legend(handles=hLst, loc="lower right", ncol=1)
        if title:
            fig.suptitle(title)
        _maybe_savefig(fig, outpath)

    return fig, axs


def plot_secfunc(
    data: dict,
    *,
    title: Optional[str] = None,
    ax: Optional[Axes] = None,
    outpath: Optional[PathLike] = None,
):
    """
    Plot the real and imaginary parts of secular functions.

    :param    data:               One dictionary from the list returned by
                                  :func:`read_secfunc <pygrt.utils.read_secfunc>`
                                  (index the list before passing), or a mapping
                                  from panel title to such dictionaries
                                  (drawn side by side).
    :param    title:              Panel title when ``data`` is a single dictionary.
    :param    ax:                 Existing matplotlib Axes. Only valid for a single dataset.
    :param    outpath:            If given, save the figure to this path.

    :return:
            - **fig** -                        matplotlib.Figure object
            - **axs** -                        tuple of matplotlib.Axes
    """
    def plot_on_ax(data: dict, ax: Axes, title: str):
        carr = np.asarray(data["c"])
        ax.plot(carr, np.asarray(data["real"]), "k-", lw=1, marker="o", ms=0, label="Real")
        ax.plot(carr, np.asarray(data["imag"]), "k--", lw=1, marker="o", ms=0, label="Imaginary")
        ax.set_xmargin(0)
        ax.set_xlabel("Phase velocity (km/s)")
        ax.set_ylabel("Secular function")
        if title:
            ax.set_title(title)
        ax.legend(loc="upper right", handletextpad=0.3, borderaxespad=0.3, labelspacing=0.2)

    if isinstance(data, (list, tuple)):
        raise TypeError("plot_secfunc() expects a dictionary; index the list returned by read_secfunc().")

    panels = _as_modal_panels(data, title, required=("c", "real", "imag"))
    fig, axs = _prepare_modal_axes(panels, ax, figsize=(5.0 * max(len(panels), 1), 3.0))
    for (panel_title, panel_data), axi in zip(panels, axs):
        plot_on_ax(panel_data, axi, panel_title)
    _maybe_savefig(fig, outpath)
    return fig, axs
