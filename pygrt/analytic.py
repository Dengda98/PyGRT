"""
    :file:     analytic.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2024-07-24

    Lamb 问题与 Okada 静态形变的解析解

"""

from __future__ import annotations

import os
import warnings
from enum import IntFlag
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
import numpy.ctypeslib as npct

from .cli import _dynamic_geometry_options, _resolve_rcv_points, format_float, format_range, run_grt


PathLike = Union[str, os.PathLike]

_NPCT_REAL_TYPE = "f8"


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

    An sCoulomb-format finite-fault file can be passed through ``src_fault``.
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
    :param    rcv_fault:        sCoulomb-format finite receiver-fault file with Kode=100. Each
                                fault contributes its center, or subfault centers
                                when ``rcv_fault_size`` is supplied. Slip magnitude is ignored;
                                an exact ``rake`` header preserves the angle even at zero slip,
                                otherwise the slip columns define direction
    :param    rcv_fault_size:   Optional ``(dL, dW)`` receiver subdivision size
                                in km along strike / dip. Both values must be positive or zero.
                                Omit or set ``(0, 0)`` to use each fault center
    :param    output_path:      Output NetCDF file path
    :param    scale:            Point-source scale in dyne·cm unless
                                ``scale_with_mu`` is true. Not used for finite faults
    :param    scale_with_mu:    If true, pass ``-Su`` and treat ``scale`` as potency
                                or area times slip in cm^3
    :param    strike:           Fault strike in degrees, in [0, 360]
    :param    dip:              Fault dip in degrees, in [0, 90]
    :param    rake:             Slip rake in degrees, in [-180, 180]
    :param    src_fault:        sCoulomb-format finite-fault file with 11 data columns;
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
    from .c_interfaces import C_grt_lamb_parse_phase_list

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

    from .c_interfaces import C_grt_solve_lamb1

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
    u = np.zeros((nt, 3, 3), dtype=_NPCT_REAL_TYPE)

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
        tbar = tbar.astype(_NPCT_REAL_TYPE, copy=False)
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

    from .c_interfaces import C_grt_solve_lamb2

    nu, tbar, R, depsrc, deprcv, azimuth = _prepare_lamb2_inputs(nu, tbar, R, depsrc, deprcv, azimuth)
    supported = _LambPhase.P | _LambPhase.S | (_LambPhase.SP if depsrc > 0 else _LambPhase.PS)
    phase_mask = _prepare_lamb_phase_mask(phases, supported)
    nt = len(tbar)
    G = np.zeros((nt, 3, 3), dtype=_NPCT_REAL_TYPE)
    Gr = np.zeros((nt, 3, 3, 3), dtype=_NPCT_REAL_TYPE)
    Gs = np.zeros((nt, 3, 3, 3), dtype=_NPCT_REAL_TYPE)
    Grs = np.zeros((nt, 3, 3, 3, 3), dtype=_NPCT_REAL_TYPE)

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

    from .c_interfaces import C_grt_solve_lamb3

    nu, tbar, R, depsrc, deprcv, azimuth = _prepare_lamb3_inputs(nu, tbar, R, depsrc, deprcv, azimuth)
    supported = _LambPhase.P | _LambPhase.S | _LambPhase.PP | _LambPhase.SS | _LambPhase.PS | _LambPhase.SP | _LambPhase.SPS
    phase_mask = _prepare_lamb_phase_mask(phases, supported)
    nt = len(tbar)
    G = np.zeros((nt, 3, 3), dtype=_NPCT_REAL_TYPE)
    Gr = np.zeros((nt, 3, 3, 3), dtype=_NPCT_REAL_TYPE)
    Gs = np.zeros((nt, 3, 3, 3), dtype=_NPCT_REAL_TYPE)
    Grs = np.zeros((nt, 3, 3, 3, 3), dtype=_NPCT_REAL_TYPE)
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
    :param    src_fault:        sCoulomb source file supplying source positions, mechanisms and
                                slip/potency. Append a complete ``-D`` option to each row
                                to specify its rupture process.
    :param    src_fault_size:   Along-strike/dip subdivision sizes (dL, dW) in km,
                                required for rectangular sources. Both values must be positive
                                or zero; ``(0, 0)`` disables further subdivision.
                                Point Kode records remain single points.
    :param    rcv_points:       ASCII receiver file: north east depth in km, optionally followed
                                by strike dip rake in degrees.
    :param    rcv_fault:        sCoulomb receiver file; only Kode=100 is supported and
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
