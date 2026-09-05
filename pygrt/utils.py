"""
    :file:     utils.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包含一些数据处理操作上的补充:   

"""

from __future__ import annotations

import os
import glob
from copy import deepcopy
from pathlib import Path
from typing import List, Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from obspy import Stream, read
from scipy.interpolate import interpn
from scipy.io import netcdf_file
from scipy.signal import oaconvolve
from scipy.special import jv
import numpy.ctypeslib as npct

from .cli import format_float, format_range, run_grt
from .c_interfaces import (
    C_grt_solve_lamb1,
    C_grt_solve_lamb2,
)


__all__ = [
    "read_static_nc",
    "read_static_grn",
    "okada",
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
    "stream_convolve",
    "stream_integral",
    "stream_diff",
    "stream_write_sac",
    "read_kernels_freqs",
    "read_statsfile",
    "read_statsfile_ptam",
    "plot_statsdata",
    "plot_statsdata_ptam",
    "lamb1",
    "solve_lamb1",
    "lamb2",
]


PathLike = Union[str, os.PathLike]


QWV_NUM = 3
INTEG_NUM = 4
SRC_M_NUM = 6
SRC_M_ORDERS = [0, 0, 1, 0, 1, 2]
SRC_M_NAME_ABBR = ["EX", "VF", "HF", "DD", "DS", "SS"]
qwvchs = ["q", "w", "v"]
NPCT_REAL_TYPE = "f8"
NPCT_CMPLX_TYPE = "c16"


def _attribute_value(value):
    """Convert a NetCDF attribute to a convenient Python value."""
    if isinstance(value, np.ndarray) and value.ndim == 0:
        return value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def read_static_nc(path: PathLike) -> dict:
    """
    Read a static NetCDF grid produced by ``grt static`` modules.

    The returned dictionary contains three top-level entries:

    * ``dimensions`` - mapping from dimension name to length
    * ``variables`` - mapping from variable name to
      ``{"dimensions", "data", "attributes"}``
    * ``attributes`` - global NetCDF attributes

    Variable arrays are available at ``variables[name]["data"]``.

    :param    path:               Path to the static NetCDF file.

    :return: A dictionary containing the NetCDF data and metadata.
    """
    path = str(path)
    if not Path(path).is_file():
        raise FileNotFoundError(f"NetCDF file does not exist: {path}")

    with netcdf_file(path, mode="r", mmap=False) as dataset:
        dimensions = {name: int(length) for name, length in dataset.dimensions.items()}
        attributes = {name: _attribute_value(getattr(dataset, name)) for name in dataset._attributes}
        variables = {}
        result = {
            "dimensions": dimensions,
            "variables": variables,
            "attributes": attributes,
        }
        for name, variable in dataset.variables.items():
            data = np.array(variable[:], copy=True)
            variable_attributes = {key: _attribute_value(value) for key, value in variable._attributes.items()}
            variables[name] = {
                "dimensions": tuple(variable.dimensions),
                "data": data,
                "attributes": variable_attributes,
            }
    return result


def read_static_grn(path: PathLike) -> dict:
    """
    Read a static Green's function NetCDF file.

    This is an alias of :func:`read_static_nc`.

    :param    path:               Path to the static Green's function file.

    :return: A dictionary containing the NetCDF data and metadata.
    """
    return read_static_nc(path)


def okada(
    *,
    modelparams: Sequence[float],
    depsrc: Optional[float] = None,
    deprcv: Optional[float] = None,
    norths: Optional[Sequence[float]] = None,
    easts: Optional[Sequence[float]] = None,
    recv_points: Optional[PathLike] = None,
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
    return_result: bool = False,
):
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
                               when ``recv_points`` is used
    :param    norths:           North grid range ``(start, stop, step)`` in km
    :param    easts:            East grid range ``(start, stop, step)`` in km
    :param    recv_points:      ASCII receiver file with either ``north east depth``
                                or ``north east depth strike dip rake``; coordinates
                                are in km and angles are in degrees
    :param    rcv_fault:        Coulomb-format finite receiver-fault file. Each
                                fault contributes its center, or subfault centers
                                when ``rcv_fault_size`` is supplied
    :param    rcv_fault_size:   Optional ``(dL, dW)`` receiver subdivision size
                                in km along strike / dip
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
    :param    return_result:    If true, read and return the generated NetCDF data

    :return: The result from :func:`read_static_nc` when ``return_result`` is true;
             otherwise ``None``
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
    vp, vs, rho = modelparams

    use_ff = src_fault is not None
    use_q = recv_points is not None
    use_r = rcv_fault is not None
    use_xy = norths is not None or easts is not None
    has_strike = strike is not None
    has_dip = dip is not None
    has_rake = rake is not None
    has_geometry = has_strike or has_dip or has_rake
    has_point_source_options = scale is not None or depsrc is not None or has_geometry or scale_with_mu

    def source_option() -> Optional[str]:
        if not has_geometry:
            return None
        if not has_strike or not has_dip:
            raise ValueError("strike and dip must be supplied together.")
        if has_rake:
            return f"-M{format_float(strike)}/{format_float(dip)}/{format_float(rake)}"
        return f"-M{format_float(strike)}/{format_float(dip)}"

    if ((use_q or use_r) and use_xy):
        raise ValueError("recv_points/rcv_fault is mutually exclusive with norths/easts.")
    if (use_q and use_r):
        raise ValueError("recv_points and rcv_fault are mutually exclusive.")
    if use_xy and (norths is None or easts is None):
        raise ValueError("norths and easts must be supplied together.")
    if ((not use_q) and (not use_r) and (not use_xy)):
        raise ValueError("Specify norths/easts, recv_points or rcv_fault.")
    if depsrc is not None and depsrc < 0.0:
        raise ValueError("depsrc must be nonnegative.")
    if deprcv is not None and deprcv < 0.0:
        raise ValueError("deprcv must be nonnegative.")
    if ((use_q or use_r) and (deprcv is not None)):
        raise ValueError("recv_points/rcv_fault is mutually exclusive with deprcv.")
    if rcv_fault_size is not None:
        if (not use_r):
            raise ValueError("rcv_fault_size requires rcv_fault.")
        if ((len(rcv_fault_size) != 2)
                or (rcv_fault_size[0] <= 0.0)
                or (rcv_fault_size[1] <= 0.0)):
            raise ValueError("rcv_fault_size must contain positive (dL, dW).")

    command = [
        "okada",
        f"-I{format_float(vp)}/{format_float(vs)}/{format_float(rho)}",
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
    elif (not use_q) and (not use_r):
        raise ValueError("deprcv is required for grid receivers.")

    if use_q:
        command.append(f"-Q{Path(recv_points)}")
    elif use_r:
        receiver_option = f"-R{Path(rcv_fault)}"
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
    if return_result:
        return read_static_nc(output)
    return None


def compute_okada(*args, **kwargs):
    """Legacy interface renamed to :func:`okada`; calling it raises an error."""
    raise RuntimeError("compute_okada() has been renamed to okada(); use okada() instead.")


def _run_static_file_module(
    path: PathLike,
    module: str,
    options: Sequence[object],
    return_result: bool,
):
    """运行只处理静态 NetCDF 文件的模块，并按需读取处理结果"""
    path = Path(path)
    if path.is_dir():
        raise ValueError(f"Only static synthesis files are supported: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Synthesis result does not exist: {path}")

    run_grt([module, f"-G{path}", *options])
    return read_static_nc(path) if return_result else None


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
        raise ValueError("lat0 and lon0 must be finite numbers.") from None
    if not np.isfinite(lat0) or not (-90.0 < lat0 < 90.0):
        raise ValueError("lat0 must be finite and in (-90, 90).")
    if not np.isfinite(lon0) or not (-180.0 <= lon0 <= 180.0):
        raise ValueError("lon0 must be finite and in [-180, 180].")

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
    recv_points: Optional[PathLike] = None,
    force_rake: bool = False,
    return_result: bool = False,
):
    """
    Project static stress tensors onto receiver-fault geometry in place.

    The input ``path`` must be a static synthesis NetCDF file, not a dynamic
    synthesis directory, and must contain the
    six stress components produced by ``static_stress``. For grid and ordinary
    points layouts, pass ``strike``, ``dip`` and ``rake`` together. For finite
    receiver points, pass only ``rake`` when the file has undefined rake values;
    set ``force_rake=True`` to replace every rake. ``recv_points`` corresponds
    to the C module's ``-Q`` option and must contain six columns per row.

    Results are written back to ``path`` as ``sigma_n`` and ``tau_s``.

    :param    path:         Static synthesis NetCDF file containing stress components.
    :param    strike:       Manual receiver strike in degrees.
    :param    dip:          Manual receiver dip in degrees.
    :param    rake:         Manual receiver rake in degrees.
    :param    recv_points:  Six-column receiver geometry file for the ``-Q`` option.
    :param    force_rake:   If true, append ``+f`` and force the manual rake for all finite points.
    :param    return_result: If true, return the processed NetCDF data.

    :return: The result from :func:`read_static_nc` when ``return_result`` is true;
             otherwise ``None``
    """
    options = []
    geometry = [value for value in (strike, dip, rake) if value is not None]
    if geometry:
        geometry_text = "/".join(format_float(value) for value in geometry)
        options.append(f"-M{geometry_text}{'+f' if force_rake else ''}")

    if recv_points is not None:
        options.append(f"-Q{Path(recv_points)}")

    return _run_static_file_module(path, "static_sproj", options, return_result)


def compute_sproj(*args, **kwargs):
    """Legacy interface renamed to :func:`static_sproj`; calling it raises an error."""
    raise RuntimeError("compute_sproj() has been renamed to static_sproj(); use static_sproj() instead.")


def static_coulomb(
    path: PathLike,
    friction: float,
    *,
    return_result: bool = False,
):
    """
    Compute Coulomb stress change in a static synthesis NetCDF file.

    The input file must already contain ``sigma_n`` and ``tau_s``, normally
    produced by :func:`static_sproj`. The result ``coulomb`` is written back
    to the same file using ``tau_s + friction * sigma_n``.

    :param    path:          Static synthesis NetCDF file containing ``sigma_n`` and ``tau_s``.
    :param    friction:      Nonnegative dimensionless effective friction coefficient.
    :param    return_result: If true, return the processed NetCDF data.

    :return: The result from :func:`read_static_nc` when ``return_result`` is true;
             otherwise ``None``
    """
    return _run_static_file_module(path, "static_coulomb", [f"-F{format_float(friction)}"], return_result)


def compute_coulomb(*args, **kwargs):
    """Legacy interface renamed to :func:`static_coulomb`; calling it raises an error."""
    raise RuntimeError("compute_coulomb() has been renamed to static_coulomb(); use static_coulomb() instead.")


def _run_dynamic_file_module(path: PathLike, module: str, return_result: bool):
    """Run a dynamic SAC-directory module and optionally read its output"""
    path = Path(path)
    if path.is_file():
        raise ValueError(f"Only dynamic synthesis directories are supported: {path}")
    if not path.is_dir():
        raise FileNotFoundError(f"Synthesis result does not exist: {path}")

    run_grt([module, path])
    return read(str(path / f"{module}_*.sac")) if return_result else None


def _run_static_tensor_module(path: PathLike, module: str, return_result: bool):
    """Run a static tensor module and optionally read its NetCDF output"""
    path = Path(path)
    if path.is_dir():
        raise ValueError(f"Only static synthesis files are supported: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Synthesis result does not exist: {path}")

    run_grt([module, path])
    return read_static_nc(path) if return_result else None


def strain(
    path: PathLike,
    *,
    return_result: bool = False,
):
    """
    Compute a dynamic strain tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same SAC directory.

    :param    path:               Dynamic SAC synthesis directory.
    :param    return_result:      If true, read and return ``strain_*.sac``.

    :return: An :class:`obspy.Stream` when ``return_result`` is true;
             otherwise ``None``.
    """
    return _run_dynamic_file_module(path, "strain", return_result)


def static_strain(
    path: PathLike,
    *,
    return_result: bool = False,
):
    """
    Compute a static strain tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same NetCDF file.

    :param    path:               Static synthesis NetCDF file.
    :param    return_result:      If true, read and return the processed NetCDF data.

    :return: A NetCDF dictionary when ``return_result`` is true; otherwise ``None``.
    """
    return _run_static_tensor_module(path, "static_strain", return_result)


def rotation(
    path: PathLike,
    *,
    return_result: bool = False,
):
    """
    Compute a dynamic rotation tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same SAC directory.

    :param    path:               Dynamic SAC synthesis directory.
    :param    return_result:      If true, read and return ``rotation_*.sac``.

    :return: An :class:`obspy.Stream` when ``return_result`` is true;
             otherwise ``None``.
    """
    return _run_dynamic_file_module(path, "rotation", return_result)


def static_rotation(
    path: PathLike,
    *,
    return_result: bool = False,
):
    """
    Compute a static rotation tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same NetCDF file.

    :param    path:               Static synthesis NetCDF file.
    :param    return_result:      If true, read and return the processed NetCDF data.

    :return: A NetCDF dictionary when ``return_result`` is true; otherwise ``None``.
    """
    return _run_static_tensor_module(path, "static_rotation", return_result)


def stress(
    path: PathLike,
    *,
    return_result: bool = False,
):
    """
    Compute a dynamic stress tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same SAC directory. Stress unit is dyne/cm² (= 0.1 Pa).

    :param    path:               Dynamic SAC synthesis directory.
    :param    return_result:      If true, read and return ``stress_*.sac``.

    :return: An :class:`obspy.Stream` when ``return_result`` is true;
             otherwise ``None``.
    """
    return _run_dynamic_file_module(path, "stress", return_result)


def static_stress(
    path: PathLike,
    *,
    return_result: bool = False,
):
    """
    Compute a static stress tensor in place from synthetic spatial derivatives.

    The synthesis must have been computed with ``calc_upar=True``. Results are
    written back to the same NetCDF file. Stress unit is dyne/cm² (= 0.1 Pa).

    :param    path:               Static synthesis NetCDF file.
    :param    return_result:      If true, read and return the processed NetCDF data.

    :return: A NetCDF dictionary when ``return_result`` is true; otherwise ``None``.
    """
    return _run_static_tensor_module(path, "static_stress", return_result)


def compute_strain(*args, **kwargs):
    """Legacy interface split into :func:`strain` and :func:`static_strain`; calling it raises an error."""
    raise RuntimeError("compute_strain() was replaced by strain() or static_strain(); use the matching interface instead.")


def compute_rotation(*args, **kwargs):
    """Legacy interface split into :func:`rotation` and :func:`static_rotation`; calling it raises an error."""
    raise RuntimeError("compute_rotation() was replaced by rotation() or static_rotation(); use the matching interface instead.")


def compute_stress(*args, **kwargs):
    """Legacy interface split into :func:`stress` and :func:`static_stress`; calling it raises an error."""
    raise RuntimeError("compute_stress() was replaced by stress() or static_stress(); use the matching interface instead.")


def stream_convolve(st0: Stream, signal0: np.ndarray, inplace: bool = True) -> Stream:
    """
    Convolve every trace with a discrete signal.

    :param    st0:            Input ObsPy stream.
    :param    signal0:        Discrete convolution signal.
    :param    inplace:        Whether to modify ``st0`` in place.

    :return: The convolved ObsPy stream.
    """
    st = st0 if inplace else deepcopy(st0)
    signal = np.asarray(signal0, dtype=float)
    for trace in st:
        dt = trace.stats.delta
        data = trace.data
        if hasattr(trace.stats, "sac") and "user0" in trace.stats.sac:
            npts = trace.stats.npts
            w_i = trace.stats.sac["user0"]
            factor = np.exp(np.arange(npts) * dt * w_i)
            adjusted_signal = signal / factor[: len(signal)]
            data[:] /= factor
            data1 = np.pad(data, (len(signal) - 1, 0), mode="wrap")
            data[:] = oaconvolve(data1, adjusted_signal, mode="valid")[:npts] * dt
            data[:] *= factor
        else:
            data1 = np.pad(data, (len(signal) - 1, 0), mode="wrap")
            data[:] = oaconvolve(data1, signal, mode="valid")[: len(data)] * dt
    return st


def stream_integral(st0: Stream, inplace: bool = True) -> Stream:
    """
    Integrate every trace with the trapezoidal rule.

    :param    st0:            Input ObsPy stream.
    :param    inplace:        Whether to modify ``st0`` in place.

    :return: The integrated ObsPy stream.
    """
    st = st0 if inplace else deepcopy(st0)
    for trace in st:
        dt = trace.stats.delta
        data = trace.data
        last = data[0]
        data[0] = 0.0
        for index in range(1, len(data)):
            current = data[index]
            data[index] = 0.5 * (current + last) * dt + data[index - 1]
            last = current
    return st


def stream_diff(st0: Stream, inplace: bool = True) -> Stream:
    """
    Differentiate every trace with a centered finite difference.

    :param    st0:            Input ObsPy stream.
    :param    inplace:        Whether to modify ``st0`` in place.

    :return: The differentiated ObsPy stream.
    """
    st = st0 if inplace else deepcopy(st0)
    for trace in st:
        trace.data[:] = np.gradient(trace.data, trace.stats.delta)
    return st


def stream_write_sac(st: Stream, directory: PathLike) -> None:
    """
    Write each trace to ``directory/{channel}.sac``.

    :param    st:             ObsPy stream to write.
    :param    directory:      Directory for the SAC files.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for trace in st:
        trace.write(str(directory / f"{trace.stats.channel}.sac"), format="SAC")


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

        :param    statsfile:       File path (Wildcards can be used to simplify input)

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

    # 获得震中距
    dist = float(os.path.dirname(statsfile).split("_")[-1])

    # 从文件路径命名中，获得对应的K文件路径
    PTAMname = os.path.basename(statsfile)
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






def lamb1(*, nu: float, tbar: np.ndarray, azimuth: float):
    r"""
        solve the first-kind Lamb's problem using the generalized closed-form solution, see：

            张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社

        :param      nu:         Poisson ratio in (0, 0.5)
        :param      tbar:       dimensionless time :math:`\bar{t}=\dfrac{t}{T_S}=\dfrac{t}{r/\beta}=\dfrac{\beta t}{r}`,
                                where :math:`T_S=r/\beta` is the S-wave time scale and :math:`r` is the direct source-receiver distance
        :param      azimuth:    azimuth in degree, in ``[0, 360]``

        :return:    Dimensionless step-force Green function with shape (nt, 3, 3). To get
                    the physical solution, divide by :math:`\pi^2 \mu r`, where
                    :math:`\mu` is the shear modulus. The returned result is dimensionless.
    """

    nu = _prepare_lamb_scalar(nu, "nu")
    tbar = _prepare_lamb_time_series(tbar)
    azimuth = _prepare_lamb_scalar(azimuth, "azimuth")
    if nu <= 0.0 or nu >= 0.5:
        raise ValueError("nu should be in (0, 0.5).")
    if azimuth < 0.0 or azimuth > 360.0:
        raise ValueError("azimuth should be in [0, 360].")

    # 定义结果数组
    nt = len(tbar)
    u = np.zeros((nt, 3, 3), dtype=NPCT_REAL_TYPE)

    C_grt_solve_lamb1(nu, npct.as_ctypes(tbar), nt, azimuth, npct.as_ctypes(u.ravel()))

    return u


def solve_lamb1(*args, **kwargs):
    """Legacy interface renamed to :func:`lamb1`; calling it raises an error."""
    raise RuntimeError("solve_lamb1() has been renamed to lamb1(); use lamb1() instead.")


def _prepare_lamb_scalar(value, name):
    try:
        value = np.asarray(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} should be a finite real scalar.") from None
    if value.ndim != 0:
        raise ValueError(f"{name} should be a finite real scalar.")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{name} should be a finite real scalar.") from None
    if not np.isfinite(value):
        raise ValueError(f"{name} should be finite.")
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
    if not np.all(np.isfinite(tbar)):
        raise ValueError("tbar should contain only finite values.")
    if np.any(tbar < 0.0):
        raise ValueError("tbar should be nonnegative.")
    if tbar.size > 1 and np.any(np.diff(tbar) <= 0.0):
        raise ValueError("tbar should be strictly increasing.")
    return np.ascontiguousarray(tbar)


def _prepare_lamb2_inputs(ts, theta, azimuth):
    ts = np.asarray(ts)
    if np.any(ts < 0.0):
        raise ValueError("ts should be nonnegative.")
    if theta <= 0.0 or theta >= 90.0:
        raise ValueError("theta should be in (0, 90) degree.")
    if azimuth < 0.0 or azimuth > 360.0:
        raise ValueError("azimuth should be in [0, 360].")
    return ts.astype(NPCT_REAL_TYPE)


def lamb2(nu:float, ts:np.ndarray, theta:float, azimuth:float):
    r"""
        Solve the second-kind Lamb problem using the generalized closed-form solution.

        The receiver is on the free surface and the source is underground. See:

            张海明, 冯禧 著. 2024. 地震学中的 Lamb 问题（下）. 科学出版社

        :param      nu:       Poisson ratio in ``(0, 0.5)``
        :param      ts:       dimensionless time :math:`\bar{t}`
        :param      theta:    source ray angle in degree, measured from the upward vertical
        :param      azimuth:  azimuth in degree, from source to receiver
        :return:    A tuple ``(G, dG_source, dG_receiver)``. ``G`` has shape
                    ``(nt, 3, 3)`` and both derivative arrays have shape
                    ``(nt, 3, 3, 3)`` with index order
                    ``[time, coordinate, receiver_component, source_component]``.
                    The normalized arrays satisfy ``G = pi^2 * mu * r * G^H``,
                    ``dG_source = pi^2 * mu * r^2 * G^H_(,k')`` and
                    ``dG_receiver = pi^2 * mu * r^2 * G^H_(,k)``.
    """

    ts = _prepare_lamb2_inputs(ts, theta, azimuth)
    nt = len(ts)
    G = np.zeros((nt, 3, 3), dtype=NPCT_REAL_TYPE)
    dG_source = np.zeros((nt, 3, 3, 3), dtype=NPCT_REAL_TYPE)
    dG_receiver = np.zeros((nt, 3, 3, 3), dtype=NPCT_REAL_TYPE)

    C_grt_solve_lamb2(
        nu,
        npct.as_ctypes(ts),
        nt,
        theta,
        azimuth,
        npct.as_ctypes(G.ravel()),
        npct.as_ctypes(dG_source.ravel()),
        npct.as_ctypes(dG_receiver.ravel()),
    )

    return G, dG_source, dG_receiver
