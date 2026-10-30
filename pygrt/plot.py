"""
    :file:     plot.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2024-07-24

    计算结果与断层几何的辅助绘图及配色

"""

from __future__ import annotations

import os
from typing import Optional, Sequence, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from scipy.special import jv

from .io import _SRC_M_NAME_ABBR, _SRC_M_ORDERS, read_sCoulomb_faults


PathLike = Union[str, os.PathLike]


def _get_stats_Fname(statsdata:np.ndarray, karr:np.ndarray, dist:float, srctype:str, ptype:str):
    # 根据ptype获得对应的核函数
    krarr = karr*dist

    # 从数组中找到震源名称的索引
    try:
        _idx = _SRC_M_NAME_ABBR.index(srctype)
        mtype = str(_SRC_M_ORDERS[_idx])
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
        Based on the data read by the :func:`read_statsfile <pygrt.io.read_statsfile>` function,
        plot the kernel function :math:`F(k,\omega)`, the integrand :math:`F(k,\omega)J_m(kr)k`,
        and calculate the cumulative integral :math:`\sum F(k,\omega)J_m(kr)k` .

        .. note:: Not every source type corresponds to every order and every integration type, see :ref:`grn_types` for details.

        :param    statsdata:         return value of :func:`read_statsfile <pygrt.io.read_statsfile>` function
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
        Based on data read by the :func:`read_statsfile_ptam <pygrt.io.read_statsfile_ptam>` function,
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

    :param    data:               Dictionary from :func:`read_dispersion <pygrt.io.read_dispersion>`,
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
                                  :func:`read_eigenfunction <pygrt.io.read_eigenfunction>`,
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
                                  :func:`read_sensitivity <pygrt.io.read_sensitivity>`.
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
                                  :func:`read_secfunc <pygrt.io.read_secfunc>`
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


def plot_sCoulomb_faults_geometry(input_path: PathLike, *, outpath: Optional[PathLike] = None):
    """
    Plot horizontal projection and 3D geometry of sCoulomb-format faults.

    Coordinates follow sCoulomb conventions: X east, Y north and depth downward, in km.
    Fault records are read using :func:`~pygrt.io.read_sCoulomb_faults`;
    row-end text is ignored when plotting.

    :param    input_path:  sCoulomb-format fault file with two header lines
    :param    outpath:     If given, save the figure to this path
    :return:              Figure and a tuple of horizontal and 3D Axes
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from matplotlib.ticker import MaxNLocator

    def fault_geometry(fault: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """计算断层顶边、底边及其三维顶点"""

        east_begin, north_begin = fault["east_begin"], fault["north_begin"]
        east_end, north_end = fault["east_end"], fault["north_end"]
        dip, top, bot = fault["dip"], fault["top"], fault["bot"]
        if dip <= 0.0 or dip > 90.0:
            raise ValueError("fault dip must be in (0, 90]")
        if bot <= top:
            raise ValueError("fault bot must be greater than top")
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

    faults = read_sCoulomb_faults(input_path)
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
        top, bot = fault["top"], fault["bot"]
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

    depth_max = max(fault["bot"] for fault in faults)
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
    ax_3d.yaxis.set_major_locator(MaxNLocator(nbins=3))
    ax_3d.set_aspect("equal")
    ax_3d.view_init(elev=24.0, azim=-58.0)
    ax_3d.legend(loc="upper left", fontsize=8)

    _maybe_savefig(fig, outpath)
    return fig, (ax_map, ax_3d)


def get_anatolia_cmap():
    """
    Return the Coulomb Anatolia colormap, from blue through white to red.

    Colors follow the Coulomb_anatolia.cpt palette used in the gallery.
    The midpoint is white, suitable for symmetric ranges around zero.

    :return:  A Matplotlib LinearSegmentedColormap
    """
    from matplotlib.colors import LinearSegmentedColormap

    # 内嵌 CPT 的 RGB 节点，使色标不依赖外部文件
    colors = np.array([
        (54, 39, 253), (55, 48, 253), (55, 57, 253),
        (56, 66, 253), (56, 74, 253), (57, 83, 253),
        (57, 92, 253), (58, 101, 253), (59, 110, 253),
        (59, 119, 253), (60, 128, 253), (60, 137, 253),
        (61, 145, 254), (61, 154, 254), (62, 163, 254),
        (62, 172, 254), (63, 181, 254), (64, 190, 254),
        (64, 199, 254), (65, 208, 254), (65, 216, 254),
        (66, 225, 254), (66, 234, 254), (67, 243, 254),
        (91, 245, 254), (114, 246, 254), (137, 247, 254),
        (161, 249, 255), (184, 251, 255), (208, 252, 255),
        (231, 253, 255), (255, 255, 255), (255, 253, 233),
        (255, 251, 211), (254, 249, 189), (254, 247, 167),
        (254, 245, 145), (253, 243, 123), (253, 241, 101),
        (253, 239, 79), (252, 229, 76), (252, 218, 74),
        (251, 208, 71), (250, 197, 69), (250, 187, 66),
        (249, 177, 64), (248, 166, 61), (248, 156, 59),
        (247, 145, 56), (246, 135, 54), (246, 125, 51),
        (245, 114, 49), (245, 104, 46), (244, 94, 44),
        (243, 83, 41), (243, 73, 39), (242, 62, 36),
        (241, 52, 34), (241, 42, 31), (240, 31, 29),
        (239, 21, 26), (239, 10, 24), (238, 0, 21),
    ], dtype=float) / 255.0
    return LinearSegmentedColormap.from_list("anatolia", colors, N=257)


def plot_static_displacement(
    data: dict, *, src_fault: Optional[PathLike] = None, title: Optional[str] = None,
    clim_ratio: float = 0.5, outpath: Optional[PathLike] = None,
):
    """
    Quickly inspect static ZNE displacement on a horizontal grid or at arbitrary points.

    Color shows vertical displacement and arrows show horizontal displacement.
    This is a quick auxiliary plot using the Anatolia colormap.
    Color limits are ±clim_ratio times the maximum absolute vertical displacement.
    Values outside this symmetric range saturate at the endpoint colors.
    For nonzero horizontal displacement, a reference arrow above the axes uses
    Matplotlib's automatic quiver scale, with a length of 1 cm or a multiple of 10 cm.

    :param    data:       Variable dictionary from :func:`~pygrt.io.read_nc_variables`, with ZNE components
    :param    src_fault:  Optional sCoulomb fault file whose top edges are drawn in black
    :param    title:      Optional plot title
    :param    clim_ratio: Color-range fraction of the maximum absolute vertical displacement, in (0, 1]
    :param    outpath:    If given, save the figure to this path
    :return:             Figure and Axes
    """
    north, east = data["north"], data["east"]
    z, n, e = data["Z"], data["N"], data["E"]
    if not np.isfinite(clim_ratio) or not 0 < clim_ratio <= 1:
        raise ValueError("clim_ratio must be finite and in (0, 1].")
    maximum = max(float(np.max(np.abs(z))) * clim_ratio, np.finfo(float).eps)
    vmin, vmax = -maximum, maximum
    cmap = get_anatolia_cmap()
    below, above = np.min(z) < vmin, np.max(z) > vmax
    extend = ("both" if below else "max") if above else ("min" if below else "neither")
    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)

    # 网格与任意点共用坐标约定，超出色标范围的值饱和显示
    if np.ndim(z) == 2:
        mesh = ax.pcolormesh(east, north, z, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
    elif np.ndim(z) == 1:
        mesh = ax.scatter(east, north, c=z, s=180, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0.0)
    else:
        raise ValueError("Displacement must be a horizontal grid or a one-dimensional point array.")
    quiver = None
    if np.any(np.hypot(e, n) > 0):
        quiver = ax.quiver(east, north, e, n, color="black", angles="xy", scale_units="xy", scale=None, width=0.003)
    if src_fault is not None:
        edges = np.loadtxt(src_fault, skiprows=2, usecols=(1, 2, 3, 4), ndmin=2)
        for x0, y0, x1, y1 in edges:
            ax.plot([x0, x1], [y0, y1], color="black", lw=3)
    ax.set_xlabel("East (km)")
    ax.set_ylabel("North (km)")
    ax.set_aspect("equal")
    ax.grid(lw=0.4)
    if title:
        ax.set_title(title, pad=38)
    fig.colorbar(mesh, ax=ax, label="Vertical displacement Z (cm)", shrink=0.7, aspect=25, extend=extend)

    if quiver is not None:
        # 首次绘制后才能取得自动 scale，参考箭头以约 15% 的图幅宽度为目标
        fig.canvas.draw()
        target = 0.15 * abs(np.diff(ax.get_xlim())[0]) * quiver.scale
        reference = max(1.0, 10.0 * np.floor(target / 10.0 + 0.5))
        ax.quiverkey(quiver, 0.78, 1.04, reference, "", coordinates="axes", labelpos="E")
        # 单独的文字参与布局与裁切，避免图外的参考标注被 tight 保存方式截掉
        ax.text(0.80, 1.04, f"{reference:g} cm", transform=ax.transAxes, va="center")

    _maybe_savefig(fig, outpath)
    return fig, ax


def plot_static_coulomb(
    data: dict, *, src_fault: Optional[PathLike] = None, clim_ratio: float = 0.5, outpath: Optional[PathLike] = None,
):
    """
    Quickly inspect shear, normal and Coulomb stress changes on a horizontal grid, in MPa.

    This is a quick auxiliary plot using a shared Anatolia color scale.
    Color limits are ±clim_ratio times the maximum absolute stress across all three components.
    Values outside this symmetric range saturate at the endpoint colors.

    :param    data:       Variable dictionary from :func:`~pygrt.io.read_nc_variables` after static Coulomb calculation
    :param    src_fault:  Optional sCoulomb fault file whose top edges are drawn in black
    :param    clim_ratio: Shared color-range fraction of the maximum absolute stress, in (0, 1]
    :param    outpath:    If given, save the figure to this path
    :return:             Figure and a tuple of three Axes
    """
    names = ("tau_s", "sigma_n", "coulomb")
    titles = ("Shear Stress Change", "Normal Stress Change", "Coulomb Stress Change")
    values = [np.asarray(data[name]) * 1e-7 for name in names]
    if not np.isfinite(clim_ratio) or not 0 < clim_ratio <= 1:
        raise ValueError("clim_ratio must be finite and in (0, 1].")
    maximum = max(max(float(np.max(np.abs(value))) for value in values) * clim_ratio, np.finfo(float).eps)
    vmin, vmax = -maximum, maximum
    cmap = get_anatolia_cmap()
    below = any(np.min(value) < vmin for value in values)
    above = any(np.max(value) > vmax for value in values)
    extend = ("both" if below else "max") if above else ("min" if below else "neither")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout='constrained')
    edges = None if src_fault is None else np.loadtxt(src_fault, skiprows=2, usecols=(1, 2, 3, 4), ndmin=2)

    # 将 dyn/cm² 转为 MPa，三个分量使用共同的色标范围
    for ax, value, title in zip(axes, values, titles):
        mesh = ax.pcolormesh(data["east"], data["north"], value, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        if edges is not None:
            for x0, y0, x1, y1 in edges:
                ax.plot([x0, x1], [y0, y1], color="black", lw=3)
        ax.set_title(title)
        ax.set_xlabel("East (km)")
        ax.set_ylabel("North (km)")
        ax.set_aspect("equal")
    fig.colorbar(mesh, ax=axes, orientation="horizontal", label="MPa", shrink=0.5, aspect=25, extend=extend)
    _maybe_savefig(fig, outpath)
    return fig, tuple(axes)
