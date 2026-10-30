"""
    :file:     cli.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2026-8-11

    调用 grt 命令行程序的工具

"""

from __future__ import annotations

import os
import warnings
from contextlib import contextmanager
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from typing import Iterable, List, Optional, Sequence, Union

import numpy as np


def find_grt() -> str:
    """Return the nearest available ``grt`` executable path."""
    bundled = Path(__file__).resolve().parent / "C_extension" / "bin" / "grt"
    if bundled.is_file():
        return str(bundled)

    from_path = shutil.which("grt")
    if from_path:
        return from_path

    raise FileNotFoundError(
        "The grt executable was not found. Build it with `make -C pygrt/C_extension CC=gcc-16 -j`."
    )


def run_grt(
    args: Sequence[object],
    *,
    cwd: Optional[Union[str, os.PathLike]] = None,
    print_log: bool = True,
) -> None:
    """
    Run one ``grt`` module and raise a detailed error on failure.

    :param    args:      Arguments passed to the ``grt`` command.
    :param    cwd:       Working directory used to run the command.
    :param    print_log: If true, stream ``grt`` stdout/stderr to the terminal.
                         If false, suppress regular stdout logs but still print
                         warnings and diagnostics. On failure, both output
                         streams are printed and included in the raised error.
    """
    command = [find_grt(), *(str(arg) for arg in args)]
    if print_log:
        completed = subprocess.run(command, cwd=cwd, check=False, stdout=sys.stdout, stderr=sys.stderr)
        if completed.returncode != 0:
            raise RuntimeError(f"grt command failed with exit code {completed.returncode}: {' '.join(command)}")
        return

    completed = subprocess.run(command, cwd=cwd, check=False, text=True, encoding="utf-8", capture_output=True)
    stdout = completed.stdout or ""
    stderr = completed.stderr or ""
    if completed.returncode != 0:
        _write_output(stdout, sys.stdout)
        _write_output(stderr, sys.stderr)
        message = f"grt command failed with exit code {completed.returncode}: {' '.join(command)}"
        detail = _format_output(stdout, stderr)
        if detail:
            message += f"\n{detail}"
        raise RuntimeError(message)

    _write_warnings(stdout)
    _write_output(stderr, sys.stderr)


def _write_output(output: str, stream) -> None:
    """Write captured process output without changing its formatting."""
    if output:
        stream.write(output)
        stream.flush()


def _write_warnings(stdout: str) -> None:
    """Forward warning lines from captured ``grt`` stdout."""
    for line in stdout.splitlines(keepends=True):
        if "[WARNING]" in line:
            _write_output(line, sys.stdout)


def _format_output(stdout: str, stderr: str) -> str:
    """Format captured process output for a Python exception."""
    details = []
    if stdout.strip():
        details.append(f"stdout:\n{stdout.rstrip()}")
    if stderr.strip():
        details.append(f"stderr:\n{stderr.rstrip()}")
    return "\n".join(details)


def format_float(value: float) -> str:
    """
    Format a numerical CLI value without unnecessary trailing zeros.

    :param    value: Numerical value to format.

    :return: Formatted command-line value.
    """
    return format(float(value), ".15g")


def format_range(values: Iterable[float], name: str) -> str:
    """
    Format three values for a ``-X`` or ``-Y`` option.

    :param    values: Three values corresponding to one CLI coordinate option.
    :param    name:   Name used in the validation error message.

    :return: Slash-separated command-line value.
    """
    try:
        values = list(values)
    except TypeError:
        raise ValueError(f"{name} must contain exactly three values.")
    if len(values) != 3:
        raise ValueError(f"{name} must contain exactly three values.")

    return "/".join(format_float(value) for value in values)


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


def _resolve_rcv_points(rcv_points: Optional[Union[str, os.PathLike]], kwargs: dict, function_name: str):
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
