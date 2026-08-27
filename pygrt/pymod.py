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
FloatOrSequence = Union[float, Sequence[float]]

__all__ = ["PyModel1D"]


def _normalize_float_array(
    values: FloatOrSequence,
    name: str,
) -> np.ndarray:
    """
    将标量或一维浮点序列规范为 float64 数组

    标量输入保持为 0 维数组，序列输入保持为 1 维数组
    """
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must be a float or a 1-D sequence of floats, not a string.")
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim == 0:
        arr = np.asarray(arr, dtype=np.float64)
    elif arr.ndim == 1:
        arr = np.ascontiguousarray(arr, dtype=np.float64)
    else:
        raise ValueError(f"{name} must be a scalar or a 1-D sequence of floats.")

    if arr.size == 0:
        raise ValueError(f"{name} must not be empty.")
    if np.any(arr < 0.0):
        raise ValueError(f"{name} must be nonnegative.")
    if arr.size > 1 and not np.all(np.diff(arr) > 0.0):
        raise ValueError(f"{name} must be strictly ascending when multiple values are given.")
    return arr


def _format_depth_list(depths: np.ndarray) -> str:
    """将深度数组格式化为 CLI ``-Ds``/``-Dr`` 的逗号列表"""
    return ",".join(format_float(float(z)) for z in depths)


class PyModel1D:
    """
    File-based 1D layered model for GRT calculations.

    Typical workflow:

    1. Create :class:`PyModel1D` with the Green's function path(s) actually needed
       (``grn`` and/or ``stgrn``) and optional ``modelpath``.
    2. Compute Green's functions with :meth:`greenfn` or :meth:`static_greenfn`
       (requires ``modelpath``).
    3. Synthesize waveforms or static fields with :meth:`syn` or
       :meth:`static_syn` (only the corresponding GF path is required).
    """

    def __init__(
        self,
        *,
        grn: Optional[PathLike] = None,
        stgrn: Optional[PathLike] = None,
        modelpath: Optional[PathLike] = None,
        topbound: str = "free",
        botbound: str = "halfspace",
    ):
        """
        Create a file-based 1D layered model handle.

        The model file is a plain text table. Each row is one layer in the form
        ``thickness(km)  Vp(km/s)  Vs(km/s)  Rho(g/cm^3)  [Qp  Qs]``.
        A zero thickness marks a half-space bottom layer.

        All arguments must be passed by keyword. ``grn`` / ``stgrn`` / ``modelpath``
        are optional at construction; methods that need them raise if missing.

        :param    grn:                Root directory for dynamic Green's functions.
        :param    stgrn:              NetCDF file path for static Green's functions.
        :param    modelpath:          Path to the layered model file. Required when
                                      computing Green's functions or travel times
                                      (unless ``modelpath`` is passed to
                                      :meth:`travt`).
        :param    topbound:           Top boundary condition. One of ``free``, ``rigid`` and ``halfspace``.
        :param    botbound:           Bottom boundary condition. One of ``free``, ``rigid`` and ``halfspace``.
        """
        if topbound not in {"free", "rigid", "halfspace"}:
            raise ValueError(f"Unsupported topbound={topbound}.")
        if botbound not in {"free", "rigid", "halfspace"}:
            raise ValueError(f"Unsupported botbound={botbound}.")

        self.topbound = topbound
        self.botbound = botbound
        self.modelpath: Optional[str] = None
        self.grn: Optional[str] = None
        self.stgrn: Optional[str] = None

        if modelpath is not None:
            self.modelpath = str(Path(modelpath))
            if not Path(self.modelpath).is_file():
                raise FileNotFoundError(f"Model file does not exist: {self.modelpath}")

        if grn is not None:
            target = Path(grn)
            target.mkdir(parents=True, exist_ok=True)
            self.grn = str(target)

        if stgrn is not None:
            target = Path(stgrn)
            target.parent.mkdir(parents=True, exist_ok=True)
            self.stgrn = str(target)

    def travt(
        self,
        *,
        depsrc: float,
        deprcv: float,
        dists: FloatOrSequence,
        modelpath: Optional[PathLike] = None,
    ):
        r"""
        Compute first-arrival P- and S-wave travel times.

        Calls the C routine ``grt_compute_travt1d_from_file``, which reads the
        layered model and evaluates travel times at the given source/receiver
        depths and epicentral distances. All arguments must be passed by keyword.

        :param    depsrc:            Source depth in km.
        :param    deprcv:            Receiver depth in km.
        :param    dists:             Epicentral distance(s) in km. A scalar or a
                                     strictly ascending sequence of distances.
        :param    modelpath:         Model file for this call only. If omitted, uses
                                     ``self.modelpath``. If both are set and differ,
                                     a warning is issued and ``self.modelpath`` is not changed.
                                     Required when ``self.modelpath`` is unset.

        :return: ``(travtP, travtS)`` in s. For a scalar distance both are floats;
                 for multiple distances both are NumPy arrays of shape ``(n,)``.
        """
        if depsrc < 0 or deprcv < 0:
            raise ValueError("Source and receiver depths must be nonnegative.")

        if modelpath is not None:
            use_model = str(Path(modelpath))
            if not Path(use_model).is_file():
                raise FileNotFoundError(f"Model file does not exist: {use_model}")
            if self.modelpath is not None:
                warnings.warn(
                    f"travt uses temporary modelpath={use_model!r}; "
                    f"instance modelpath={self.modelpath!r} is unchanged.",
                    UserWarning,
                    stacklevel=2,
                )
        else:
            if self.modelpath is None:
                raise RuntimeError("Pass modelpath= to travt() or set it in PyModel1D(...).")
            use_model = self.modelpath

        # 标量震中距返回 float，序列返回长度为 n 的数组
        distances = _normalize_float_array(dists, "dists")
        single = (distances.ndim == 0)
        distances = np.atleast_1d(distances)

        carr = C_grt_compute_travt1d_from_file(
            use_model.encode("utf-8"),
            float(depsrc),
            float(deprcv),
            distances.ctypes.data_as(PREAL),
            c_size_t(distances.size),
        )
        if cast(carr, c_void_p).value is None:
            raise RuntimeError(f"Failed to compute travel times for model {use_model}.")

        arr = npct.as_array(carr, shape=(distances.size, 2)).copy()
        C_grt_free(carr)

        if single:
            return float(arr[0, 0]), float(arr[0, 1])
        return arr[:, 0].copy(), arr[:, 1].copy()

    def compute_travt1d(self, *args, **kwargs):
        """Legacy interface renamed to :meth:`travt`; calling it raises an error."""
        raise RuntimeError("compute_travt1d() has been renamed to travt(); use travt() instead.")

    def greenfn(
        self,
        *,
        depsrc: FloatOrSequence,
        deprcv: FloatOrSequence,
        dists: FloatOrSequence,
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

        Requires ``grn`` and ``modelpath`` in the constructor. Results are written
        as SAC files under ``{grn}/{model}_{depsrc}_{deprcv}_{distance}``. The
        model file is also copied into the ``grn`` root directory by the C module.
        All arguments must be passed by keyword.

        :param    depsrc:            Source depth or strictly ascending source-depth
                                     sequence in km.
        :param    deprcv:            Receiver depth or strictly ascending
                                     receiver-depth sequence in km.
        :param    dists:             Array of strictly ascending epicentral distances
                                     in km, or a single distance.
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
        if self.grn is None:
            raise RuntimeError("Pass grn= to PyModel1D(...) before greenfn().")
        if self.modelpath is None:
            raise RuntimeError("Pass modelpath= to PyModel1D(...) before greenfn().")
        if nt <= 0 or dt <= 0:
            raise ValueError("nt and dt must be positive.")
        depsrcs = np.atleast_1d(_normalize_float_array(depsrc, "depsrc"))
        deprcvs = np.atleast_1d(_normalize_float_array(deprcv, "deprcv"))
        distances = np.atleast_1d(_normalize_float_array(dists, "dists"))

        try:
            freq1, freq2 = freqband
        except (TypeError, ValueError):
            raise ValueError("freqband must contain exactly two values (f1, f2).") from None

        multi_depth = (depsrcs.size > 1) or (deprcvs.size > 1)
        command = {
            "subcommand": "greenfn",
            "M": f"-M{self.modelpath}",
        }
        if multi_depth:
            command["Ds"] = f"-Ds{_format_depth_list(depsrcs)}"
            command["Dr"] = f"-Dr{_format_depth_list(deprcvs)}"
        else:
            command["D"] = f"-D{format_float(float(depsrcs[0]))}/{format_float(float(deprcvs[0]))}"
        command.update({
            "N": f"-N{nt}/{format_float(dt)}+w{format_float(zeta)}+n{upsampling_n}",
            "R": f"-R{','.join(format_float(distance) for distance in distances)}",
            "O": f"-O{self.grn}",
            "B": f"-B{self._boundary_option()}",
        })

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

    def compute_grn(self, *args, **kwargs):
        """Legacy interface renamed to :meth:`greenfn`; calling it raises an error."""
        raise RuntimeError("compute_grn() has been renamed to greenfn(); use greenfn() instead.")

    def static_greenfn(
        self,
        *,
        depsrc: FloatOrSequence,
        deprcv: FloatOrSequence,
        dists: Optional[FloatOrSequence] = None,
        norths: Optional[Sequence[float]] = None,
        easts: Optional[Sequence[float]] = None,
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
        Compute static Green's functions with the ``grt static_greenfn`` command.

        Requires ``stgrn`` and ``modelpath`` in the constructor. Results are written
        to the configured NetCDF file (4D STGRNLIB layout
        ``[depsrc][deprcv][north][east]``) and currently overwrite any existing
        content. All arguments must be passed by keyword.

        ``depsrc`` / ``deprcv`` may be a scalar or a 1-D sequence:

        * Single depth pair: CLI ``-Ddepsrc/deprcv``.
        * Multiple depths: CLI ``-Ds...`` / ``-Dr...`` (comma-separated list).

        For an ordinary Green's-function library, ``dists`` is the usual
        choice: it computes a one-dimensional distance list and stores it as
        ``north=0`` and ``east=dists``. This is also the most convenient
        format for :meth:`static_syn`, which synthesizes results at
        surrounding distance samples and combines those results with
        distance-based weights.
        Receiver locations can be specified in either of two ways:

        1. ``dists``, a scalar or sequence of epicentral distances in km.
           This is equivalent to placing receivers along the east axis with
           north = 0 and is the recommended way to build a reusable library.
        2. ``norths`` and ``easts``, each a three-value sequence
           ``(start, stop, step)`` in km, mapped to CLI ``-X`` / ``-Y``. This
           form is mainly useful for accuracy tests or when a later synthesis
           should reuse the same grid without distance-based combination of
           multiple synthesized results.

        :param    depsrc:            Source depth(s) in km. Multiple values must be
                                     strictly ascending.
        :param    deprcv:            Receiver depth(s) in km. Multiple values must be
                                     strictly ascending.
        :param    dists:             Epicentral distance(s) in km. A scalar or
                                     strictly ascending sequence, equivalent to
                                     receivers with north = 0 and east = ``dists``.
                                     Mutually exclusive with ``norths`` / ``easts``.
        :param    norths:            Three values defining the north-coordinate
                                     option ``-Xstart/stop/step`` in km.
        :param    easts:             Three values defining the east-coordinate
                                     option ``-Ystart/stop/step`` in km.
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
        if self.stgrn is None:
            raise RuntimeError("Pass stgrn= to PyModel1D(...) before static_greenfn().")
        if self.modelpath is None:
            raise RuntimeError("Pass modelpath= to PyModel1D(...) before static_greenfn().")

        depsrcs = np.atleast_1d(_normalize_float_array(depsrc, "depsrc"))
        deprcvs = np.atleast_1d(_normalize_float_array(deprcv, "deprcv"))
        multi_depth = (depsrcs.size > 1) or (deprcvs.size > 1)

        if dists is not None:
            if norths is not None or easts is not None:
                raise ValueError("Use either dists or norths/easts.")
            distances = np.atleast_1d(_normalize_float_array(dists, "dists"))
            command_grid = {
                "R": f"-R{','.join(format_float(value) for value in distances)}"
            }
        else:
            if norths is None or easts is None:
                raise ValueError("Set norths and easts, or set dists.")
            command_grid = {
                "X": f"-X{format_range(norths, 'norths')}",
                "Y": f"-Y{format_range(easts, 'easts')}",
            }

        command = {
            "module": "static_greenfn",
            "M": f"-M{self.modelpath}",
        }
        if multi_depth:
            command["Ds"] = f"-Ds{_format_depth_list(depsrcs)}"
            command["Dr"] = f"-Dr{_format_depth_list(deprcvs)}"
        else:
            command["D"] = f"-D{format_float(float(depsrcs[0]))}/{format_float(float(deprcvs[0]))}"
        command["O"] = f"-O{self.stgrn}"
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

    def compute_static_grn(self, *args, **kwargs):
        """Legacy interface renamed to :meth:`static_greenfn`; calling it raises an error."""
        raise RuntimeError("compute_static_grn() has been renamed to static_greenfn(); use static_greenfn() instead.")

    def syn(
        self,
        *,
        depsrc: Optional[float] = None,
        deprcv: Optional[float] = None,
        dist: Optional[float] = None,
        azimuth: float,
        output_path: PathLike,
        scale: float,
        scale_with_mu: bool = False,
        strike: Optional[float] = None,
        dip: Optional[float] = None,
        rake: Optional[float] = None,
        force: Optional[Sequence[float]] = None,
        moment_tensor: Optional[Sequence[float]] = None,
        time_function: Optional[str] = None,
        integrate_order: Optional[int] = None,
        differentiate_order: Optional[int] = None,
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

        Call :meth:`greenfn` first (or point ``grn`` at an existing GF
        root or subdirectory). When ``grn`` is a root, a selector is required
        for a dimension with multiple values. For a singleton dimension it may
        be omitted or explicitly set, but an explicit value must match. When
        ``grn`` is a subdirectory, all three selectors must be omitted. No
        interpolation is performed. All arguments must be passed by keyword.

        The source type is inferred from the source-specific parameters. Leave
        ``strike``, ``dip``, ``rake``, ``force`` and ``moment_tensor`` unset for
        an explosion (``EX``). Supplying ``force`` selects a single force
        (``SF``); supplying ``moment_tensor`` selects a moment tensor (``MT``);
        supplying ``strike`` and ``dip`` selects a tensile crack (``TS``), or a
        double-couple (``DC``) when ``rake`` is also supplied. Only one source
        parameter group may be used at a time.

        :param    depsrc:              Source depth in km when ``grn`` is a GF root
                                       with multiple source depths. It may be omitted
                                       when the root has one source depth, but an
                                       explicit value must match the root library.
                                       It must be omitted when ``grn`` is a subdirectory.
        :param    deprcv:              Receiver depth in km when ``grn`` is a GF root
                                       with multiple receiver depths. It may be omitted
                                       when the root has one receiver depth, but an
                                       explicit value must match the root library.
                                       It must be omitted when ``grn`` is a subdirectory.
        :param    dist:                Epicentral distance in km when ``grn`` is a GF
                                       root with multiple distances. It may be omitted
                                       when the root has one distance, but an explicit
                                       value must match the root library. It must be
                                       omitted when ``grn`` is a subdirectory.
        :param    azimuth:             Azimuth from source to receiver in deg.
                                       North is 0°, clockwise positive.
        :param    output_path:         Output directory for SAC files
                                       ``{output_path}/{ch}.sac``.
        :param    scale:               Source scaling factor. For explosion,
                                       double-couple, tensile-crack and moment-tensor
                                       sources, this is the scalar seismic moment in
                                       dyne·cm. For a single force, the unit is dyne.
                                       If ``scale_with_mu`` is true, ``scale`` is
                                       treated as area × slip in cm³ and multiplied by
                                       the source-layer shear modulus :math:`\mu`.
        :param    scale_with_mu:       If true, multiply ``scale`` by the source-layer
                                       shear modulus :math:`\mu` (CLI ``-Su``).
        :param    strike:              Fault strike in deg, in [0, 360]. North is 0°,
                                       clockwise positive. Set together with ``dip``
                                       for ``TS`` or with ``dip`` and ``rake`` for
                                       ``DC``.
        :param    dip:                 Fault dip in deg, in [0, 90]. Required for
                                       ``TS`` and ``DC`` when the corresponding
                                       source geometry is selected.
        :param    rake:                Slip rake in deg, in [-180, 180],
                                       counterclockwise positive on the fault plane.
                                       Supplying it with ``strike`` and ``dip`` selects
                                       ``DC``; omit it for ``TS``.
        :param    force:               Single-force coefficients ``(fN, fE, fZ)`` for
                                       the ``SF`` source. Directions are north, east
                                       and downward. Each coefficient is multiplied by
                                       ``scale``.
        :param    moment_tensor:       Six independent moment-tensor coefficients
                                       ``(Mxx, Mxy, Mxz, Myy, Myz, Mzz)`` for the
                                       ``MT`` source. Subscripts x/y/z denote
                                       north/east/down.
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
        if self.grn is None:
            raise RuntimeError("Pass grn= to PyModel1D(...) before syn().")
        if dist is not None and dist < 0:
            raise ValueError("dist must be nonnegative.")
        output = Path(output_path)
        output.mkdir(parents=True, exist_ok=True)

        command = {
            "subcommand": "syn",
            "G": f"-G{self.grn}",
            "A": f"-A{format_float(azimuth)}",
            "S": f"-S{'u' if scale_with_mu else ''}{format_float(scale)}",
            "O": f"-O{output}",
        }
        if depsrc is not None:
            if depsrc < 0:
                raise ValueError("depsrc must be nonnegative.")
            command["Ds"] = f"-Ds{format_float(depsrc)}"
        if deprcv is not None:
            if deprcv < 0:
                raise ValueError("deprcv must be nonnegative.")
            command["Dr"] = f"-Dr{format_float(deprcv)}"
        if dist is not None:
            command["R"] = f"-R{format_float(dist)}"
        command.update(self._source_options(strike, dip, rake, force, moment_tensor))

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

    def compute_syn(self, *args, **kwargs):
        """Legacy interface renamed to :meth:`syn`; calling it raises an error."""
        raise RuntimeError("compute_syn() has been renamed to syn(); use syn() instead.")

    def static_syn(
        self,
        *,
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
        force: Optional[Sequence[float]] = None,
        moment_tensor: Optional[Sequence[float]] = None,
        src_fault: Optional[PathLike] = None,
        src_fault_size: Optional[Sequence[float]] = None,
        zne: bool = False,
        calc_upar: bool = False,
        return_result: bool = False,
    ):
        r"""
        Synthesize static three-component displacement with ``grt static_syn``.

        Results are written to the NetCDF file ``output_path``. Requires
        ``stgrn`` (and typically a prior :meth:`static_greenfn`).
        All arguments must be passed by keyword.

        Receivers default to the library north/east grid. Optionally redefine
        them with ``norths``/``easts`` (uniform ``deprcv`` when the library has
        multiple receiver depths), or with ``recv_points`` for an ASCII file of
        arbitrary ``north east depth`` points (CLI ``-Q``). Each row may append
        ``strike dip rake`` in degrees; these angles are saved in the output but
        are not used in synthesis. ``recv_points`` is
        mutually exclusive with ``norths``/``easts`` and ``deprcv``. If the
        library was built with ``dists`` / ``-R``, the default grid is a 1-D
        line (north = 0, east = R); set ``norths``/``easts`` or ``recv_points``
        to obtain a 2-D field.
        A Coulomb-format finite receiver-fault file can be supplied through
        ``rcv_fault`` (CLI ``-R``); without ``rcv_fault_size`` the subdivision
        size defaults to the smallest positive sampling interval among epicentral
        distance, source depth and receiver depth in the library. With
        ``rcv_fault_size``, each fault is subdivided with ``(dL, dW)`` along
        strike / dip.

        Point-source type is inferred from the source-specific parameters:
        leaving ``strike``, ``dip``, ``rake``, ``force`` and ``moment_tensor``
        unset selects an explosion (``EX``); ``force`` selects a single force
        (``SF``); ``moment_tensor`` selects a moment tensor (``MT``); ``strike``
        and ``dip`` select a tensile crack (``TS``), or a double-couple (``DC``)
        when ``rake`` is also supplied. Only one source parameter group may be
        used at a time. Finite faults use ``src_fault`` (Coulomb-format file,
        CLI ``-C``) instead. Its ``Kode`` column selects rectangular shear/
        tensile sources or point shear/expansion sources. An exact ``rake``
        token in the seventh header column selects Kode 100 rake/net-slip
        rows; the filename suffix is not used to select the format. That path
        requires a multi-source-depth library and rejects point-source options.

        For each target receiver, the C module first synthesizes results at the
        surrounding epicentral-distance samples and combines those synthesized
        results with weights based on the target distance. When a requested
        source or receiver depth lies between samples, it performs the same
        process for each surrounding depth combination and then combines those
        synthesized results with depth-based weights. Thus, interpolation is
        applied to synthesized results rather than directly to Green's-function
        arrays. If the library was generated with an explicit ``-X``/``-Y`` grid
        and the same grid is reused, the corresponding synthesized results can
        be used directly.

        :param    depsrc:            Point-source depth in km (CLI ``-Ds``). Required
                                     when the library has multiple source depths;
                                     optional when it has one, but an explicit value
                                     must match the library. Forbidden when
                                     ``src_fault`` is set.
        :param    deprcv:            Receiver depth in km for grid receivers
                                     (CLI ``-Dr``). Required when the library has
                                     multiple receiver depths and ``recv_points`` is
                                     not used; optional when it has one, but an
                                     explicit value must match the library. Do not set
                                     it when using ``recv_points``.
        :param    norths:            Optional new north grid as three values
                                     ``(start, stop, step)`` in km. Must be set
                                     together with ``easts``. Mutually exclusive
                                     with ``recv_points``.
        :param    easts:             Optional new east grid as three values
                                     ``(start, stop, step)`` in km. Must be set
                                     together with ``norths``. Mutually exclusive
                                     with ``recv_points``.
        :param    recv_points:       ASCII file of arbitrary receivers
                                     (``north east depth`` in km, optionally
                                     followed by ``strike dip rake`` in degrees;
                                     ``#`` comments). All data rows must use the
                                     same 3- or 6-column format. Mutually exclusive
                                     with ``norths``/``easts`` and ``deprcv``.
        :param    rcv_fault:        Coulomb-format finite receiver-fault file with 11 data
                                     columns; an exact ``rake`` token in the seventh header
                                     column selects Kode 100 rake/net-slip interpretation
                                     (CLI ``-R``). Without ``rcv_fault_size``,
                                     the library sampling intervals determine the
                                     default subdivision size. With that argument,
                                     each fault contributes multiple subfault
                                     centers. Mutually exclusive with
                                     ``recv_points``, ``norths``/``easts`` and
                                     ``deprcv``.
        :param    rcv_fault_size:   Optional positive ``(dL, dW)`` in km for
                                     receiver-fault subdivision along strike / dip;
                                     if omitted, use the smallest positive interval
                                     among epicentral distance, source depth and
                                     receiver depth in the library.
        :param    output_path:       Output NetCDF file path.
        :param    scale:             Point-source scaling factor. For explosion,
                                     double-couple, tensile-crack and moment-tensor
                                     sources, this is the scalar seismic moment in
                                     dyne·cm. For a single force, the unit is dyne.
                                     If ``scale_with_mu`` is true, ``scale`` is
                                     treated as area × slip in cm³ and multiplied by
                                     the source-layer shear modulus :math:`\mu`.
                                     Required for point sources; ignored for
                                     ``src_fault``.
        :param    scale_with_mu:     If true, multiply ``scale`` by the source-layer
                                     shear modulus :math:`\mu` (CLI ``-Su``).
        :param    strike:            Fault strike in deg, in [0, 360]. North is 0°,
                                     clockwise positive. Set together with ``dip``
                                     for ``TS`` or with ``dip`` and ``rake`` for
                                     ``DC``.
        :param    dip:               Fault dip in deg, in [0, 90]. Required for
                                     ``TS`` and ``DC`` when the corresponding
                                     source geometry is selected.
        :param    rake:              Slip rake in deg, in [-180, 180],
                                     counterclockwise positive on the fault plane.
                                     Supplying it with ``strike`` and ``dip`` selects
                                     ``DC``; omit it for ``TS``.
        :param    force:             Single-force coefficients ``(fN, fE, fZ)`` for
                                     the ``SF`` source. Directions are north, east
                                     and downward. Each coefficient is multiplied by
                                     ``scale``.
        :param    moment_tensor:     Six independent moment-tensor coefficients
                                     ``(Mxx, Mxy, Mxz, Myy, Myz, Mzz)`` for the
                                     ``MT`` source. Subscripts x/y/z denote
                                     north/east/down.
        :param    src_fault:         Coulomb-format finite-fault file with 11 data columns
                                     (Kode 100/200/300/400/500). An exact ``rake`` token in
                                     the seventh header column selects Kode 100 rake/net-slip
                                     interpretation; the filename suffix is not used. Mutually
                                     exclusive with point-source options; point-source arguments
                                     cause ``ValueError``.
        :param    src_fault_size:    Optional ``(dL, dW)`` in km for finite-fault
                                     subdivision along strike / dip. If omitted,
                                     the C code uses the smallest positive interval
                                     among epicentral distance, source depth and
                                     receiver depth in the library.
        :param    zne:               If true, output ZNE instead of ZRT components.
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
        if self.stgrn is None:
            raise RuntimeError("Pass stgrn= to PyModel1D(...) before static_syn().")
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)

        use_ff = src_fault is not None
        use_q = recv_points is not None
        use_r = rcv_fault is not None
        use_xy = norths is not None or easts is not None
        has_geometry = strike is not None or dip is not None or rake is not None
        has_force = force is not None
        has_moment_tensor = moment_tensor is not None
        has_point_source_options = (
            scale is not None
            or scale_with_mu
            or has_geometry
            or has_force
            or has_moment_tensor
        )
        if use_ff and has_point_source_options:
            raise ValueError("src_fault is mutually exclusive with point-source options.")
        if ((use_q or use_r) and use_xy):
            raise ValueError("recv_points/rcv_fault is mutually exclusive with norths/easts.")
        if (use_q and use_r):
            raise ValueError("recv_points and rcv_fault are mutually exclusive.")
        if ((use_q or use_r) and (deprcv is not None)):
            raise ValueError("recv_points/rcv_fault is mutually exclusive with deprcv.")
        if use_xy and (norths is None or easts is None):
            raise ValueError("norths and easts must be supplied together.")
        if depsrc is not None and depsrc < 0.0:
            raise ValueError("depsrc must be nonnegative.")
        if deprcv is not None and deprcv < 0.0:
            raise ValueError("deprcv must be nonnegative.")
        if use_ff and depsrc is not None:
            raise ValueError("depsrc is forbidden when src_fault is set.")
        if src_fault_size is not None:
            if not use_ff:
                raise ValueError("src_fault_size requires src_fault.")
            if len(src_fault_size) != 2:
                raise ValueError("src_fault_size must be (dL, dW).")
        if rcv_fault_size is not None:
            if (not use_r):
                raise ValueError("rcv_fault_size requires rcv_fault.")
            if ((len(rcv_fault_size) != 2)
                    or (rcv_fault_size[0] <= 0.0)
                    or (rcv_fault_size[1] <= 0.0)):
                raise ValueError("rcv_fault_size must contain positive (dL, dW).")

        command = {
            "module": "static_syn",
            "G": f"-G{self.stgrn}",
            "O": f"-O{output}",
        }

        if use_ff:
            c_opt = f"-C{Path(src_fault)}"
            if src_fault_size is not None:
                c_opt += f"+i{format_float(src_fault_size[0])}/{format_float(src_fault_size[1])}"
            command["C"] = c_opt
        else:
            if scale is None:
                raise ValueError("scale is required for point-source synthesis.")
            command["S"] = f"-S{'u' if scale_with_mu else ''}{format_float(scale)}"
            command.update(self._source_options(strike, dip, rake, force, moment_tensor))
            if depsrc is not None:
                command["Ds"] = f"-Ds{format_float(depsrc)}"

        if deprcv is not None:
            command["Dr"] = f"-Dr{format_float(deprcv)}"
        if use_q:
            command["Q"] = f"-Q{Path(recv_points)}"
        elif use_r:
            r_opt = f"-R{Path(rcv_fault)}"
            if rcv_fault_size is not None:
                r_opt += f"+i{format_float(rcv_fault_size[0])}/{format_float(rcv_fault_size[1])}"
            command["R"] = r_opt
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

    def compute_static_syn(self, *args, **kwargs):
        """Legacy interface renamed to :meth:`static_syn`; calling it raises an error."""
        raise RuntimeError("compute_static_syn() has been renamed to static_syn(); use static_syn() instead.")


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
        strike: Optional[float],
        dip: Optional[float],
        rake: Optional[float],
        force: Optional[Sequence[float]],
        moment_tensor: Optional[Sequence[float]],
    ) -> Dict[str, str]:
        has_strike = strike is not None
        has_dip = dip is not None
        has_rake = rake is not None
        has_geometry = has_strike or has_dip or has_rake
        has_force = force is not None
        has_moment_tensor = moment_tensor is not None
        if has_force:
            if has_geometry or has_moment_tensor:
                raise ValueError("force is mutually exclusive with strike/dip/rake and moment_tensor.")
            if len(force) != 3:
                raise ValueError("force must contain exactly three values: (fN, fE, fZ).")
            return {"F": "-F" + "/".join(format_float(value) for value in force)}
        if has_moment_tensor:
            if has_geometry:
                raise ValueError("moment_tensor is mutually exclusive with strike/dip/rake and force.")
            if len(moment_tensor) != 6:
                raise ValueError(
                    "moment_tensor must contain exactly six values: (Mxx, Mxy, Mxz, Myy, Myz, Mzz)."
                )
            return {"T": "-T" + "/".join(format_float(value) for value in moment_tensor)}
        if not has_geometry:
            return {}
        if not has_strike or not has_dip:
            raise ValueError("strike and dip must be supplied together.")
        if has_rake:
            return {"M": f"-M{format_float(strike)}/{format_float(dip)}/{format_float(rake)}"}
        return {"M": f"-M{format_float(strike)}/{format_float(dip)}"}
