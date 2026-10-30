"""
    :file:     io.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2024-07-24

    计算结果、积分统计与断层文件的读取和解析

"""

from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import List, Sequence, Union

import numpy as np
from scipy.interpolate import interpn
from scipy.io import netcdf_file


PathLike = Union[str, os.PathLike]

_QWV_NUM = 3
_INTEG_NUM = 4
_SRC_M_NUM = 6
_SRC_M_ORDERS = [0, 0, 1, 0, 1, 2]
_SRC_M_NAME_ABBR = ["EX", "VF", "HF", "DD", "DS", "SS"]
_qwvchs = ["q", "w", "v"]
_NPCT_REAL_TYPE = "f8"
_NPCT_CMPLX_TYPE = "c16"


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
    dtype = [('k' if basename[0] == 'K' else 'c', _NPCT_REAL_TYPE)]
    for im in range(_SRC_M_NUM):
        modr = _SRC_M_ORDERS[im]
        for c in range(_QWV_NUM):
            if modr==0 and _qwvchs[c] == 'v':
                continue

            dtype.append((f"{_SRC_M_NAME_ABBR[im]}_{_qwvchs[c]}", _NPCT_CMPLX_TYPE))


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
    for im in range(_SRC_M_NUM):
        modr = _SRC_M_ORDERS[im]
        for v in range(_INTEG_NUM):
            if modr==0 and v!=0 and v!=2:
                continue

            dtype.append((f"sum_{_SRC_M_NAME_ABBR[im]}_{v}_k", _NPCT_REAL_TYPE))
            dtype.append((f"sum_{_SRC_M_NAME_ABBR[im]}_{v}", _NPCT_CMPLX_TYPE))


    ptam_data = np.fromfile(statsfile, dtype=dtype)

    return data1, data2, ptam_data, dist


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


def read_sCoulomb_faults(input_path: PathLike) -> list[dict]:
    """
    Read an sCoulomb-format fault file into a list of dictionaries.

    Each dictionary contains the original columns named after the C FINITE_FAULT fields:
    ``east_begin``, ``north_begin``, ``east_end``, ``north_end``, ``kode``, ``value1``,
    ``value2``, ``dip``, ``top`` and ``bot``. The first placeholder column is discarded.
    ``kode`` is an integer; the other numeric fields are floats.
    X/Y correspond to east/north, coordinates and depths are in km, and dip is in degrees.
    ``value1`` and ``value2`` retain their input values, including in rake/net-slip files.

    Text after ``bot`` is stored in ``other``, with surrounding whitespace stripped and
    internal whitespace preserved; absent text gives an empty string.
    Blank lines and full-line comments after the two headers are skipped.
    No derived quantities are computed, and geometry and slip values are not validated.

    :param    input_path:  sCoulomb-format fault file with two header lines
    :return:              A list of fault dictionaries in file order
    """
    input_path = Path(input_path)
    keys = ("east_begin", "north_begin", "east_end", "north_end", "kode", "value1", "value2", "dip", "top", "bot")
    faults = []
    with input_path.open(encoding="utf-8") as file:
        header = file.readline().split()
        if len(header) < 11 or header[0] != "#" or len(file.readline().split()) < 11:
            raise ValueError(f"{input_path}: expected two sCoulomb header lines")

        # 仅拆分前 11 个数值列，保留行末文本中的空白与未知选项
        for line_number, line in enumerate(file, start=3):
            fields = line.split(maxsplit=11)
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) < 11:
                raise ValueError(f"{input_path}:{line_number}: expected 11 numeric columns")
            try:
                values = [float(field) for field in fields[:11]]
            except ValueError as exc:
                raise ValueError(f"{input_path}:{line_number}: invalid numeric value") from exc
            if not values[5].is_integer():
                raise ValueError(f"{input_path}:{line_number}: kode must be an integer")

            fault = dict(zip(keys, values[1:]))
            fault["kode"] = int(values[5])
            fault["other"] = fields[11].strip() if len(fields) > 11 else ""
            faults.append(fault)

    if not faults:
        raise ValueError(f"{input_path}: no faults found")
    return faults
