"""
    :file:     pymod.py  
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)  
    :date:     2024-07-24  

    该文件包括 Python 端使用的基于文件的模型 :class:`PyModel1D`

"""

from __future__ import annotations

import os
import warnings
from ctypes import c_size_t, cast, c_void_p
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Union

import numpy as np
import numpy.ctypeslib as npct
from obspy import read

from .cli import format_float, format_range, run_grt
from .c_interfaces import C_grt_compute_travt1d_from_file, C_grt_free, PREAL
from .utils import read_static_nc


PathLike = Union[str, os.PathLike]
DepthLike = Union[float, Sequence[float]]

__all__ = ["PyModel1D"]


def _normalize_distarr(distarr):
    """
    将 distarr 规范为一维 float64 数组

    仅接受标量浮点数或一维浮点序列；字符串等类型直接拒绝
    """
    if isinstance(distarr, (str, bytes)):
        raise TypeError("distarr must be a float or a 1-D sequence of floats, not a string.")
    arr = np.asarray(distarr, dtype=np.float64)
    if arr.ndim == 0:
        return True, np.ascontiguousarray([float(arr)], dtype=np.float64)
    if arr.ndim == 1:
        return False, np.ascontiguousarray(arr, dtype=np.float64)
    raise ValueError("distarr must be a scalar or a 1-D sequence of floats.")


def _normalize_depths(depths: DepthLike, name: str) -> np.ndarray:
    """
    将震源/接收深度规范为一维 float64 数组

    接受标量或一维浮点序列；空数组与负深度会报错
    """
    if isinstance(depths, (str, bytes)):
        raise TypeError(f"{name} must be a float or a 1-D sequence of floats, not a string.")
    arr = np.asarray(depths, dtype=np.float64)
    if arr.ndim == 0:
        arr = np.ascontiguousarray([float(arr)], dtype=np.float64)
    elif arr.ndim == 1:
        arr = np.ascontiguousarray(arr, dtype=np.float64)
    else:
        raise ValueError(f"{name} must be a scalar or a 1-D sequence of floats.")
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty.")
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be nonnegative.")
    return arr


def _format_depth_list(depths: np.ndarray) -> str:
    """将深度数组格式化为 CLI ``-Ds``/``-Dr`` 的逗号列表"""
    return ",".join(format_float(float(z)) for z in depths)


class PyModel1D:
    """
    File-based 1D layered model for GRT calculations.

    Typical workflow:

    1. Create the model from a layered-model file.
    2. Call :meth:`set_dynamic_grn_path` or :meth:`set_static_grn_path`.
    3. Compute Green's functions with :meth:`compute_grn` or :meth:`compute_static_grn`.
    4. Synthesize waveforms or static fields with :meth:`compute_syn` or :meth:`compute_static_syn`.
    """

    def __init__(
        self,
        modelpath: PathLike,
        topbound: str = "free",
        botbound: str = "halfspace",
    ):
        """
        Create a file-based 1D layered model.

        The model file is a plain text table. Each row is one layer in the form
        ``thickness(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  [Qp  Qs]``.
        A zero thickness marks a half-space bottom layer.

        :param    modelpath:          Path to the layered model file.
        :param    topbound:           Top boundary condition. One of ``free``, ``rigid`` and ``halfspace``.
        :param    botbound:           Bottom boundary condition. One of ``free``, ``rigid`` and ``halfspace``.
        """
        self.modelpath = str(Path(modelpath))
        self.topbound = topbound
        self.botbound = botbound
        self.dynamic_grn_path: Optional[str] = None
        self.static_grn_path: Optional[str] = None

        if not Path(self.modelpath).is_file():
            raise FileNotFoundError(f"Model file does not exist: {self.modelpath}")
        if topbound not in {"free", "rigid", "halfspace"}:
            raise ValueError(f"Unsupported topbound={topbound}.")
        if botbound not in {"free", "rigid", "halfspace"}:
            raise ValueError(f"Unsupported botbound={botbound}.")

    def compute_travt1d(
        self,
        *,
        depsrc: float,
        deprcv: float,
        distarr: Union[float, Sequence[float]],
    ):
        r"""
        Compute first-arrival P- and S-wave travel times.

        Calls the C routine ``grt_compute_travt1d_from_file``, which reads the
        layered model from ``modelpath`` and evaluates travel times at the given
        source/receiver depths and epicentral distances. All arguments must be
        passed by keyword.

        :param    depsrc:            Source depth in km.
        :param    deprcv:            Receiver depth in km.
        :param    distarr:           Epicentral distance(s) in km. A scalar or a sequence of distances.

        :return: ``(travtP, travtS)`` in s. For a scalar distance both are floats;
                 for multiple distances both are NumPy arrays of shape ``(n,)``.
        """
        if depsrc < 0 or deprcv < 0:
            raise ValueError("Source and receiver depths must be nonnegative.")

        # 标量震中距返回 float，序列返回长度为 n 的数组
        single, distances = _normalize_distarr(distarr)
        if distances.size == 0 or np.any(distances < 0):
            raise ValueError("distarr must contain nonnegative distances.")

        carr = C_grt_compute_travt1d_from_file(
            self.modelpath.encode("utf-8"),
            float(depsrc),
            float(deprcv),
            distances.ctypes.data_as(PREAL),
            c_size_t(distances.size),
        )
        if cast(carr, c_void_p).value is None:
            raise RuntimeError(f"Failed to compute travel times for model {self.modelpath}.")

        arr = npct.as_array(carr, shape=(distances.size, 2)).copy()
        C_grt_free(carr)

        if single:
            return float(arr[0, 0]), float(arr[0, 1])
        return arr[:, 0].copy(), arr[:, 1].copy()

    def set_dynamic_grn_path(self, path: PathLike) -> str:
        """
        Set and create the root directory for dynamic Green's functions.

        Later calls to :meth:`compute_grn` write SAC files under this directory.
        Subdirectories are named
        ``{model}_{depsrc}_{deprcv}_{distance}``.

        :param    path:               Root directory for dynamic Green's functions.

        :return: The configured dynamic Green's function directory.
        """
        target = Path(path)
        target.mkdir(parents=True, exist_ok=True)
        self.dynamic_grn_path = str(target)
        return self.dynamic_grn_path

    def set_static_grn_path(self, path: PathLike) -> str:
        """
        Set the NetCDF file path for static Green's functions.

        Later calls to :meth:`compute_static_grn` write (and currently overwrite)
        this file. Parent directories are created if needed.

        :param    path:               NetCDF file path for static Green's functions.

        :return: The configured static Green's function file path.
        """
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        self.static_grn_path = str(target)
        return self.static_grn_path

    def compute_grn(
        self,
        *,
        depsrc: float,
        deprcv: float,
        distarr: Union[float, Sequence[float]],
        nt: int,
        dt: float,
        upsampling_n: int = 1,
        freqband: Sequence[float] = (-1.0, -1.0),
        zeta: float = 0.8,
        keepAllFreq: bool = False,
        vmin_ref: float = 0.0,
        keps: float = -1.0,
        ampk: float = 2.0,
        k0: float = 50.0,
        use_kmax_ref: bool = False,
        Length: float = 0.0,
        filonLength: float = 0.0,
        safilonTol: float = 0.0,
        filonCut: float = 0.0,
        converg_method: str = "AUTO",
        delayT0: float = 0.0,
        delayV0: float = 0.0,
        ref_first_p: bool = False,
        skipImagComps: bool = False,
        calc_upar: bool = False,
        gf_source: Optional[Iterable[str]] = None,
        statsidxs: Optional[Sequence[int]] = None,
        print_log: bool = True,
    ):
        r"""
        Compute dynamic Green's functions with the ``grt greenfn`` command.

        Call :meth:`set_dynamic_grn_path` first. Results are written as SAC files
        under ``{dynamic_grn_path}/{model}_{depsrc}_{deprcv}_{distance}``.
        All arguments must be passed by keyword.

        :param    depsrc:            Source depth in km.
        :param    deprcv:            Receiver depth in km.
        :param    distarr:           Array of epicentral distances in km, or a single distance.
        :param    nt:                Number of time points. With the help of SciPy,
                                     ``nt`` no longer needs to be a power of 2.
        :param    dt:                Time interval in s.
        :param    upsampling_n:      Upsampling factor applied after inverse FFT.
        :param    freqband:          Frequency range ``(f1, f2)`` in Hz. Negative values mean
                                     that the corresponding bound is determined automatically.
        :param    zeta:              Coefficient defining the imaginary angular frequency,
                                     :math:`\tilde{\omega} = \omega - j w_I`,
                                     where :math:`w_I = \zeta\pi/T` and :math:`T=nt\,dt`.
        :param    keepAllFreq:       Whether to calculate all frequency points,
                                     regardless of how low the frequency is.
        :param    vmin_ref:          Minimum reference velocity in km/s. ``0.0`` means the
                                     minimum model velocity, limited to 0.1 km/s.
        :param    keps:              Automatic convergence condition. See Yao and Harkrider
                                     (1983) for more details. A negative value disables this condition.
        :param    ampk:              Amplification factor in the reference maximum wavenumber.
        :param    k0:                Coefficient in the reference maximum wavenumber,
                                     :math:`k_{\mathrm{max,ref}} =
                                     \sqrt{(k_0\pi/h_s)^2 + (ampk\,\omega/v_{\mathrm{min,ref}})^2}`.
                                     Here :math:`h_s=\max(|depsrc-deprcv|,0.1)`.
        :param    use_kmax_ref:      Whether to use the reference maximum wavenumber directly,
                                     without amplitude searching.
        :param    Length:            Integration step :math:`dk=2\pi/(Lr_{\max})`.
                                     ``0.0`` means the value is selected automatically.
        :param    filonLength:       Integration step of fixed-interval Filon integration at
                                     large distances, but not at zero distance. ``0.0`` disables
                                     this method. Do not set together with ``safilonTol``.
        :param    safilonTol:        Precision of self-adaptive Filon integration at large
                                     distances, but not at zero distance. ``0.0`` disables this
                                     method. Do not set together with ``filonLength``.
        :param    filonCut:          Splitting point between DWM and Filon integration,
                                     :math:`k^*=\mathrm{filonCut}/r_{\max}`.
        :param    converg_method:    Explicit convergence method. One of ``AUTO``, ``NONE``,
                                     ``DCM`` and ``PTAM``.
        :param    delayT0:           Time delay at zero distance in s.
        :param    delayV0:           Reference velocity for the time delay in km/s.
                                     Used only when ``ref_first_p`` is false.
        :param    ref_first_p:       Whether to use the first P-wave arrival as the reference
                                     for the time delay (CLI ``-Ep``).
        :param    skipImagComps:     Whether to skip the amplitude compensation from the
                                     imaginary frequency.
        :param    calc_upar:         Whether to calculate spatial derivatives of displacement.
                                     Required later if strain, stress or rotation will be computed.
        :param    gf_source:         Source types to calculate. Choose from ``EX``, ``VF``,
                                     ``HF`` and ``DC``. ``None`` means all available source types.
        :param    statsidxs:         Frequency indexes for optional statistics output.
                                     ``None`` means no statistics files. An empty list means
                                     all frequency indexes (CLI bare ``-S``).
        :param    print_log:         Whether to print calculation logs.

        :return: ``None``. Results are written to disk.
        """
        if self.dynamic_grn_path is None:
            raise RuntimeError("Call set_dynamic_grn_path() before compute_grn().")
        if nt <= 0 or dt <= 0:
            raise ValueError("nt and dt must be positive.")
        if depsrc < 0 or deprcv < 0:
            raise ValueError("Source and receiver depths must be nonnegative.")

        _, distances = _normalize_distarr(distarr)
        if distances.size == 0 or np.any(distances < 0):
            raise ValueError("distarr must contain nonnegative distances.")

        try:
            freq1, freq2 = freqband
        except (TypeError, ValueError):
            raise ValueError("freqband must contain exactly two values (f1, f2).") from None

        command = {
            "subcommand": "greenfn",
            "M": f"-M{self.modelpath}",
            "D": f"-D{format_float(depsrc)}/{format_float(deprcv)}",
            "N": f"-N{nt}/{format_float(dt)}+w{format_float(zeta)}+n{upsampling_n}",
            "R": f"-R{','.join(format_float(distance) for distance in distances)}",
            "O": f"-O{self.dynamic_grn_path}",
            "B": f"-B{self._boundary_option()}",
        }

        # Build the -N option.
        if keepAllFreq:
            command["N"] += "+a"
        if skipImagComps:
            command["N"] += "+f"

        # Build the -H option.
        command["H"] = f"-H{format_float(freq1)}/{format_float(freq2)}"

        # Build the -L option.
        option = format_float(Length)
        if filonLength:
            option += f"+l{format_float(filonLength)}"
        if safilonTol:
            option += f"+a{format_float(safilonTol)}"
        if filonCut:
            option += f"+o{format_float(filonCut)}"
        command["L"] = f"-L{option}"

        # Build the -C option.
        option = self._convergence_option(converg_method)
        if option:
            command["C"] = f"-C{option}"

        # Build the -K option.
        options = [f"+k{format_float(k0)}"]
        if use_kmax_ref:
            options.append("+f")
        options.extend([f"+s{format_float(ampk)}", f"+e{format_float(keps)}"])
        if vmin_ref:
            options.append(f"+v{format_float(vmin_ref)}")
        command["K"] = "-K" + "".join(options)

        # Build the -E option.
        if ref_first_p:
            command["E"] = f"-Ep{format_float(delayT0)}"
        else:
            command["E"] = f"-E{format_float(delayT0)}/{format_float(delayV0)}"

        # Build the -G option.
        if gf_source is not None:
            source_codes = {"EX": "e", "VF": "v", "HF": "h", "DC": "s"}
            codes = []
            for name in gf_source:
                key = str(name).upper()
                if key not in source_codes:
                    raise ValueError(f"Unsupported gf_source={name!r}. Choose from EX, VF, HF and DC.")
                codes.append(source_codes[key])
            command["G"] = "-G" + "".join(codes)

        # Build the -S option.
        if statsidxs is not None:
            command["S"] = "-S" + ",".join(str(index) for index in statsidxs)

        # Build the derivative and logging options.
        if calc_upar:
            command["e"] = "-e"
        if not print_log:
            command["s"] = "-s"

        run_grt(list(command.values()), print_log=print_log)

    def compute_static_grn(
        self,
        *,
        depsrc: DepthLike,
        deprcv: DepthLike,
        norths: Optional[Sequence[float]] = None,
        easts: Optional[Sequence[float]] = None,
        distarr: Optional[Sequence[float]] = None,
        keps: float = -1.0,
        k0: float = 50.0,
        use_kmax_ref: bool = False,
        Length: float = 15.0,
        filonLength: float = 0.0,
        safilonTol: float = 0.0,
        filonCut: float = 0.0,
        converg_method: str = "AUTO",
        calc_upar: bool = False,
        stats: bool = False,
    ):
        r"""
        Compute static Green's functions with the ``grt static greenfn`` command.

        Call :meth:`set_static_grn_path` first. Results are written to the
        configured NetCDF file (4D STGRNLIB layout
        ``[depsrc][deprcv][north][east]``) and currently overwrite any existing
        content. All arguments must be passed by keyword.

        ``depsrc`` / ``deprcv`` may be a scalar or a 1-D sequence:

        * Single depth pair: CLI ``-Ddepsrc/deprcv``.
        * Multiple depths: CLI ``-Ds...`` / ``-Dr...`` (comma-separated list).

        Receiver locations can be specified in either of two ways:

        1. ``norths`` and ``easts``, each a three-value sequence
           ``(start, stop, step)`` in km, mapped to CLI ``-X`` / ``-Y``.
        2. ``distarr``, a list of epicentral distances in km. This is equivalent
           to placing receivers along the east axis with north = 0.

        :param    depsrc:            Source depth(s) in km. Multiple values must be
                                     strictly ascending.
        :param    deprcv:            Receiver depth(s) in km. Multiple values must be
                                     strictly ascending.
        :param    norths:            Three values defining the north-coordinate
                                     option ``-Xstart/stop/step`` in km.
        :param    easts:             Three values defining the east-coordinate
                                     option ``-Ystart/stop/step`` in km.
        :param    distarr:           Epicentral distances in km. Equivalent to
                                     receivers with north = 0 and east = ``distarr``.
                                     Mutually exclusive with ``norths`` / ``easts``.
        :param    keps:              Automatic convergence condition. See Yao and
                                     Harkrider (1983) for more details. A negative
                                     value disables this condition.
        :param    k0:                Coefficient in the reference maximum wavenumber,
                                     where :math:`h_s=\max(|depsrc-deprcv|,0.1)`.
                                     The actual maximum wavenumber is searched based
                                     on the kernel amplitude.
        :param    use_kmax_ref:      Whether to use the reference maximum wavenumber
                                     directly, without amplitude searching.
        :param    Length:            Integration step :math:`dk=2\pi/(Lr_{\max})`.
                                     The default is 15.
        :param    filonLength:       Step parameter for fixed-interval Filon
                                     integration at large distances. ``0.0`` disables
                                     this method. Do not set together with
                                     ``safilonTol``.
        :param    safilonTol:        Tolerance for self-adaptive Filon integration
                                     at large distances. ``0.0`` disables this method.
                                     Do not set together with ``filonLength``.
        :param    filonCut:          Splitting point between DWM and Filon integration.
        :param    converg_method:    Explicit convergence method. One of
                                     ``AUTO``, ``NONE``, ``DCM`` and ``PTAM``.
        :param    calc_upar:         Whether to calculate spatial derivatives of
                                     displacement. Required later if strain,
                                     stress or rotation will be computed.
        :param    stats:             Whether to write integration statistics.
                                     Only available for a single source/receiver depth;
                                     ignored with a warning for multi-depth runs.

        :return: ``None``. Results are written to the configured NetCDF file.
        """
        if self.static_grn_path is None:
            raise RuntimeError("Call set_static_grn_path() before compute_static_grn().")

        depsrcs = _normalize_depths(depsrc, "depsrc")
        deprcvs = _normalize_depths(deprcv, "deprcv")
        if depsrcs.size > 1 and not np.all(np.diff(depsrcs) > 0.0):
            raise ValueError("depsrc must be strictly ascending when multiple values are given.")
        if deprcvs.size > 1 and not np.all(np.diff(deprcvs) > 0.0):
            raise ValueError("deprcv must be strictly ascending when multiple values are given.")
        multi_depth = (depsrcs.size > 1) or (deprcvs.size > 1)

        if distarr is not None:
            if norths is not None or easts is not None:
                raise ValueError("Use either distarr or norths/easts.")
            _, distances = _normalize_distarr(distarr)
            if distances.size == 0 or np.any(distances < 0.0):
                raise ValueError("distarr must contain nonnegative distances.")
            if distances.size > 1 and not np.all(np.diff(distances) > 0.0):
                raise ValueError("distarr must be strictly ascending.")
            command_grid = {
                "R": f"-R{','.join(format_float(value) for value in distances)}"
            }
        else:
            if norths is None or easts is None:
                raise ValueError("Set norths and easts, or set distarr.")
            command_grid = {
                "X": f"-X{format_range(norths, 'norths')}",
                "Y": f"-Y{format_range(easts, 'easts')}",
            }

        command = {
            "module": "static",
            "subcommand": "greenfn",
            "M": f"-M{self.modelpath}",
        }
        if multi_depth:
            command["Ds"] = f"-Ds{_format_depth_list(depsrcs)}"
            command["Dr"] = f"-Dr{_format_depth_list(deprcvs)}"
        else:
            command["D"] = f"-D{format_float(float(depsrcs[0]))}/{format_float(float(deprcvs[0]))}"
        command["O"] = f"-O{self.static_grn_path}"
        command["B"] = f"-B{self._boundary_option()}"
        command.update(command_grid)

        # Build the -L option.
        option = format_float(Length)
        if filonLength:
            option += f"+l{format_float(filonLength)}"
        if safilonTol:
            option += f"+a{format_float(safilonTol)}"
        if filonCut:
            option += f"+o{format_float(filonCut)}"
        command["L"] = f"-L{option}"

        # Build the -C option.
        option = self._convergence_option(converg_method)
        if option:
            command["C"] = f"-C{option}"

        # Build the -K option.
        options = [f"+k{format_float(k0)}"]
        if use_kmax_ref:
            options.append("+f")
        options.append(f"+e{format_float(keps)}")
        command["K"] = "-K" + "".join(options)

        # Build the statistics and derivative options.
        if stats:
            if multi_depth:
                warnings.warn(
                    "stats is ignored for multi-depth STGRNLIB computation.",
                    stacklevel=2,
                )
            else:
                command["S"] = "-S"
        if calc_upar:
            command["e"] = "-e"

        run_grt(list(command.values()))

    def compute_syn(
        self,
        *,
        dist: float,
        azimuth: float,
        scale: float,
        output_path: PathLike,
        source: str = "EX",
        strike: Optional[float] = None,
        dip: Optional[float] = None,
        rake: Optional[float] = None,
        force: Optional[Sequence[float]] = None,
        moment_tensor: Optional[Sequence[float]] = None,
        time_function: Optional[str] = None,
        integrate_order: Optional[int] = None,
        differentiate_order: Optional[int] = None,
        scale_with_mu: bool = False,
        zne: bool = False,
        calc_upar: bool = False,
        return_result: bool = False,
    ):
        r"""
        Synthesize dynamic three-component displacement with ``grt syn``.

        Results are written as SAC files under ``output_path``. By default the
        synthetics are impulse-like displacements in cm with ZRT components:

        * ``Z`` - vertical upward
        * ``R`` - radial outward
        * ``T`` - clockwise 90° from ``R``

        Call :meth:`set_dynamic_grn_path` and :meth:`compute_grn` first. The
        Green's function directory is located under ``dynamic_grn_path`` by
        matching ``dist`` in the subdirectory name. All arguments must be
        passed by keyword.

        Choose one source type with ``source``:

        * ``EX`` - explosion. Only ``scale`` is required.
        * ``DC`` - double-couple / shear. Requires ``strike``, ``dip`` and ``rake``.
        * ``TS`` - tensile crack. Requires ``strike`` and ``dip``.
        * ``SF`` - single force. Requires ``force=(fN, fE, fZ)``.
        * ``MT`` - moment tensor. Requires
          ``moment_tensor=(Mxx, Mxy, Mxz, Myy, Myz, Mzz)``.

        :param    dist:                Epicentral distance in km. Used to locate the
                                       Green's function directory under
                                       ``dynamic_grn_path``.
        :param    azimuth:             Azimuth from source to receiver in deg.
                                       North is 0°, clockwise positive.
        :param    scale:               Source scaling factor. For ``EX``, ``DC``,
                                       ``TS`` and ``MT``, this is the scalar seismic
                                       moment in dyne·cm. For ``SF``, the unit is dyne.
                                       If ``scale_with_mu`` is true, ``scale`` is
                                       treated as area × slip in cm³ and multiplied by
                                       the source-layer shear modulus :math:`\mu`.
        :param    output_path:         Output directory for SAC files
                                       ``{output_path}/{ch}.sac``.
        :param    source:              Source type. One of ``EX``, ``DC``, ``TS``,
                                       ``SF`` and ``MT``.
        :param    strike:              Fault strike in deg, in [0, 360]. North is 0°,
                                       clockwise positive. Required for ``DC`` and
                                       ``TS``.
        :param    dip:                 Fault dip in deg, in [0, 90]. Required for
                                       ``DC`` and ``TS``.
        :param    rake:                Slip rake in deg, in [-180, 180],
                                       counterclockwise positive on the fault plane.
                                       Required for ``DC``.
        :param    force:               Single-force coefficients ``(fN, fE, fZ)`` for
                                       ``SF``. Directions are north, east and downward.
                                       Each coefficient is multiplied by ``scale``.
        :param    moment_tensor:       Six independent moment-tensor coefficients
                                       ``(Mxx, Mxy, Mxz, Myy, Myz, Mzz)`` for ``MT``.
                                       Subscripts x/y/z denote north/east/down.
        :param    time_function:       Time-function string passed to CLI ``-D``.
                                       Supported forms include:

                                       * ``p/t0`` - parabola lasting ``t0`` s
                                       * ``t/t1/t2/t3`` - trapezoid with rise,
                                         plateau and fall cutoffs in s
                                       * ``r/f0`` - Ricker wavelet with dominant
                                         frequency ``f0`` in Hz
                                       * ``0/file`` - custom one-column amplitude file

                                       The peak amplitude of the time function is 1.
                                       Omit this argument for an impulse response.
        :param    integrate_order:     Number of time integrations. For example,
                                       ``1`` yields step-like displacement.
        :param    differentiate_order: Number of time differentiations. For example,
                                       ``1`` yields velocity.
        :param    scale_with_mu:       If true, multiply ``scale`` by the source-layer
                                       shear modulus :math:`\mu` (CLI ``-Su``).
        :param    zne:                 If true, output ZNE instead of ZRT components.
        :param    calc_upar:           If true, also synthesize spatial derivatives of
                                       displacement. Derivative channel names are
                                       prefixed with ``z``, ``r`` or ``t``. Set this
                                       when strain, stress or rotation will be computed
                                       later.
        :param    return_result:       If true, read the generated SAC files into an
                                       :class:`obspy.Stream`.

        :return: An ObsPy stream when ``return_result`` is true; otherwise ``None``.
        """
        grn_path = self._dynamic_grn_dir(dist)
        output = Path(output_path)
        output.mkdir(parents=True, exist_ok=True)

        command = {
            "subcommand": "syn",
            "G": f"-G{grn_path}",
            "A": f"-A{format_float(azimuth)}",
            "S": f"-S{'u' if scale_with_mu else ''}{format_float(scale)}",
            "O": f"-O{output}",
        }
        command.update(
            self._source_options(source, strike, dip, rake, force, moment_tensor)
        )

        # Build the time-function and operation-order options.
        if time_function is not None:
            command["D"] = f"-D{time_function}"
        if integrate_order is not None:
            command["I"] = f"-I{integrate_order}"
        if differentiate_order is not None:
            command["J"] = f"-J{differentiate_order}"

        # Build the component and derivative options.
        if zne:
            command["N"] = "-N"
        if calc_upar:
            command["e"] = "-e"

        run_grt(list(command.values()))
        if return_result:
            return read(str(output / "*.sac"))
        return None

    def compute_static_syn(
        self,
        *,
        output_path: PathLike,
        scale: Optional[float] = None,
        source: str = "EX",
        strike: Optional[float] = None,
        dip: Optional[float] = None,
        rake: Optional[float] = None,
        force: Optional[Sequence[float]] = None,
        moment_tensor: Optional[Sequence[float]] = None,
        scale_with_mu: bool = False,
        depsrc: Optional[float] = None,
        deprcv: Optional[float] = None,
        norths: Optional[Sequence[float]] = None,
        easts: Optional[Sequence[float]] = None,
        recv_points: Optional[PathLike] = None,
        finite_fault: Optional[PathLike] = None,
        subfault_size: Optional[Sequence[float]] = None,
        zne: bool = False,
        calc_upar: bool = False,
        return_result: bool = False,
    ):
        r"""
        Synthesize static three-component displacement with ``grt static syn``.

        Results are written to the NetCDF file ``output_path``. Call
        :meth:`set_static_grn_path` and :meth:`compute_static_grn` first.
        All arguments must be passed by keyword.

        Receivers default to the library north/east grid. Optionally redefine
        them with ``norths``/``easts`` (uniform ``deprcv`` when the library has
        multiple receiver depths), or with ``recv_points`` for an ASCII file of
        arbitrary ``north east depth`` points (CLI ``-Q``). ``recv_points`` is
        mutually exclusive with ``norths``/``easts`` and ``deprcv``. If the
        library was built with ``distarr`` / ``-R``, the default grid is a 1-D
        line (north = 0, east = R); set ``norths``/``easts`` or ``recv_points``
        to obtain a 2-D field.

        Point sources use ``scale`` and ``source``. Finite faults use
        ``finite_fault`` (Coulomb-format file, CLI ``-C``) instead; that path
        requires a multi-source-depth library, automatically writes ZNE, and
        ignores point-source options.

        Choose one point-source type with ``source``:

        * ``EX`` - explosion. Only ``scale`` is required.
        * ``DC`` - double-couple / shear. Requires ``strike``, ``dip`` and ``rake``.
        * ``TS`` - tensile crack. Requires ``strike`` and ``dip``.
        * ``SF`` - single force. Requires ``force=(fN, fE, fZ)``.
        * ``MT`` - moment tensor. Requires
          ``moment_tensor=(Mxx, Mxy, Mxz, Myy, Myz, Mzz)``.

        :param    output_path:       Output NetCDF file path.
        :param    scale:             Point-source scaling factor. For ``EX``, ``DC``,
                                     ``TS`` and ``MT``, this is the scalar seismic
                                     moment in dyne·cm. For ``SF``, the unit is dyne.
                                     If ``scale_with_mu`` is true, ``scale`` is
                                     treated as area × slip in cm³ and multiplied by
                                     the source-layer shear modulus :math:`\mu`.
                                     Required for point sources; ignored for
                                     ``finite_fault``.
        :param    source:            Point-source type. One of ``EX``, ``DC``,
                                     ``TS``, ``SF`` and ``MT``. Ignored when
                                     ``finite_fault`` is set.
        :param    strike:            Fault strike in deg, in [0, 360]. North is 0°,
                                     clockwise positive. Required for ``DC`` and
                                     ``TS``.
        :param    dip:               Fault dip in deg, in [0, 90]. Required for
                                     ``DC`` and ``TS``.
        :param    rake:              Slip rake in deg, in [-180, 180],
                                     counterclockwise positive on the fault plane.
                                     Required for ``DC``.
        :param    force:             Single-force coefficients ``(fN, fE, fZ)`` for
                                     ``SF``. Directions are north, east and downward.
                                     Each coefficient is multiplied by ``scale``.
        :param    moment_tensor:     Six independent moment-tensor coefficients
                                     ``(Mxx, Mxy, Mxz, Myy, Myz, Mzz)`` for ``MT``.
                                     Subscripts x/y/z denote north/east/down.
        :param    scale_with_mu:     If true, multiply ``scale`` by the source-layer
                                     shear modulus :math:`\mu` (CLI ``-Su``).
        :param    depsrc:            Point-source depth in km (CLI ``-Ds``). Required
                                     when the library has multiple source depths;
                                     forbidden for ``finite_fault``.
        :param    deprcv:            Receiver depth in km for grid receivers
                                     (CLI ``-Dr``). Required when the library has
                                     multiple receiver depths and ``recv_points``
                                     is not used.
        :param    norths:            Optional new north grid as three values
                                     ``(start, stop, step)`` in km. Must be set
                                     together with ``easts``. Mutually exclusive
                                     with ``recv_points``.
        :param    easts:             Optional new east grid as three values
                                     ``(start, stop, step)`` in km. Must be set
                                     together with ``norths``. Mutually exclusive
                                     with ``recv_points``.
        :param    recv_points:       ASCII file of arbitrary receivers
                                     (``north east depth`` in km; ``#`` comments).
                                     Mutually exclusive with ``norths``/``easts``
                                     and ``deprcv``.
        :param    finite_fault:      Coulomb-format finite-fault file (CLI ``-C``).
                                     Mutually exclusive with point-source options.
        :param    subfault_size:     Optional ``(dL, dW)`` in km for finite-fault
                                     subdivision along strike / dip. If omitted,
                                     the C code uses ``min(dr, dz)`` of the library.
        :param    zne:               If true, output ZNE instead of ZRT components.
                                     Finite faults always write ZNE.
        :param    calc_upar:         If true, also synthesize spatial derivatives of
                                     displacement. Derivative variable names use
                                     prefixes ``z``/``r``/``t`` (ZRT) or
                                     ``z``/``n``/``e`` (ZNE). Set this when strain,
                                     stress or rotation will be computed later.
        :param    return_result:     If true, read the generated NetCDF file with
                                     :func:`pygrt.utils.read_static_nc`.

        :return: The synthesized NetCDF data when ``return_result`` is true;
                 otherwise ``None``.
        """
        if self.static_grn_path is None:
            raise RuntimeError("Call set_static_grn_path() before compute_static_syn().")
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)

        use_ff = finite_fault is not None
        use_q = recv_points is not None
        use_xy = norths is not None or easts is not None
        if use_ff and (
            scale is not None
            or force is not None
            or moment_tensor is not None
            or strike is not None
            or dip is not None
            or rake is not None
            or scale_with_mu
            or source.upper() != "EX"
        ):
            raise ValueError("finite_fault is mutually exclusive with point-source options.")
        if use_q and use_xy:
            raise ValueError("recv_points is mutually exclusive with norths/easts.")
        if use_q and deprcv is not None:
            raise ValueError("recv_points is mutually exclusive with deprcv.")
        if use_xy and (norths is None or easts is None):
            raise ValueError("norths and easts must be supplied together.")
        if depsrc is not None and depsrc < 0.0:
            raise ValueError("depsrc must be nonnegative.")
        if deprcv is not None and deprcv < 0.0:
            raise ValueError("deprcv must be nonnegative.")
        if use_ff and depsrc is not None:
            raise ValueError("depsrc is forbidden when finite_fault is set.")
        if subfault_size is not None:
            if not use_ff:
                raise ValueError("subfault_size requires finite_fault.")
            if len(subfault_size) != 2:
                raise ValueError("subfault_size must be (dL, dW).")

        command = {
            "module": "static",
            "subcommand": "syn",
            "G": f"-G{self.static_grn_path}",
            "O": f"-O{output}",
        }

        if use_ff:
            c_opt = f"-C{Path(finite_fault)}"
            if subfault_size is not None:
                c_opt += f"+i{format_float(subfault_size[0])}/{format_float(subfault_size[1])}"
            command["C"] = c_opt
        else:
            if scale is None:
                raise ValueError("scale is required for point-source synthesis.")
            command["S"] = f"-S{'u' if scale_with_mu else ''}{format_float(scale)}"
            command.update(
                self._source_options(source, strike, dip, rake, force, moment_tensor)
            )
            if depsrc is not None:
                command["Ds"] = f"-Ds{format_float(depsrc)}"

        if deprcv is not None:
            command["Dr"] = f"-Dr{format_float(deprcv)}"
        if use_q:
            command["Q"] = f"-Q{Path(recv_points)}"
        elif use_xy:
            command["X"] = f"-X{format_range(norths, 'norths')}"
            command["Y"] = f"-Y{format_range(easts, 'easts')}"

        if zne:
            command["N"] = "-N"
        if calc_upar:
            command["e"] = "-e"

        run_grt(list(command.values()))
        if return_result:
            return read_static_nc(output)
        return None

    def _dynamic_grn_dir(self, dist: float) -> str:
        """
        在 dynamic_grn_path 下按震中距匹配格林函数子目录

        子目录命名为 ``{model}_{depsrc}_{deprcv}_{dist}``
        当前假设仅有一套震源/台站深度，故只需匹配 dist
        """
        if self.dynamic_grn_path is None:
            raise RuntimeError("Call set_dynamic_grn_path() before compute_syn().")
        root = Path(self.dynamic_grn_path)
        if not root.is_dir():
            raise FileNotFoundError(f"Dynamic Green's function root does not exist: {root}")

        suffix = f"_{format_float(dist)}"
        matches = [path for path in root.iterdir() if path.is_dir() and path.name.endswith(suffix)]
        if not matches:
            raise FileNotFoundError(f"No Green's function directory matching dist={format_float(dist)} under {root}.")
        if len(matches) > 1:
            names = ", ".join(path.name for path in sorted(matches))
            raise RuntimeError(f"Multiple Green's function directories match dist={format_float(dist)} under {root}: {names}.")
        return str(matches[0])

    def _boundary_option(self) -> str:
        return {
            "free": "f",
            "rigid": "r",
            "halfspace": "h",
        }[self.topbound] + {
            "free": "F",
            "rigid": "R",
            "halfspace": "H",
        }[self.botbound]

    @staticmethod
    def _convergence_option(value: str) -> str:
        options = {"AUTO": "", "DCM": "d", "PTAM": "p", "NONE": "n"}
        try:
            return options[value.upper()]
        except KeyError:
            raise ValueError(f"Unsupported convergence method: {value}") from None

    @staticmethod
    def _source_options(
        source: str,
        strike: Optional[float],
        dip: Optional[float],
        rake: Optional[float],
        force: Optional[Sequence[float]],
        moment_tensor: Optional[Sequence[float]],
    ) -> Dict[str, str]:
        source = source.upper()
        if source == "EX":
            return {}
        if source == "DC":
            if strike is None or dip is None or rake is None:
                raise ValueError("DC source requires strike, dip and rake.")
            return {"M": f"-M{format_float(strike)}/{format_float(dip)}/{format_float(rake)}"}
        if source == "TS":
            if strike is None or dip is None:
                raise ValueError("TS source requires strike and dip.")
            return {"M": f"-M{format_float(strike)}/{format_float(dip)}"}
        if source == "SF":
            if force is None or len(force) != 3:
                raise ValueError("SF source requires force=(fN, fE, fZ).")
            return {"F": "-F" + "/".join(format_float(value) for value in force)}
        if source == "MT":
            if moment_tensor is None or len(moment_tensor) != 6:
                raise ValueError("MT source requires six moment-tensor values.")
            return {"T": "-T" + "/".join(format_float(value) for value in moment_tensor)}
        raise ValueError(f"Unsupported source type: {source}")
