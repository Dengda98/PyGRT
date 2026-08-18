"""
比较 C CLI 与 Python 文件工作流结果的辅助函数
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

import numpy as np
from obspy import read
from scipy.io import netcdf_file


PathLike = Union[str, Path]


def _rel_error(a: np.ndarray, b: np.ndarray) -> Optional[float]:
    """相对误差；双方全零时返回 None 表示跳过"""
    if np.all(a == 0.0) and np.all(b == 0.0):
        return None
    denom = np.mean(np.abs(b))
    if denom == 0.0:
        denom = np.mean(np.abs(a))
    if denom == 0.0:
        return None
    return float(np.sum(np.abs(a - b)) / denom)


def compare_sac_dirs(
    dir_py: PathLike,
    dir_c: PathLike,
    pattern: str = "*.sac",
) -> float:
    """
    比较两个目录下同名 SAC 波形，返回平均相对误差

    以文件名（而非 channel）配对。应变/旋转/应力的 channel 可能同名
    （如 EE），只能靠 strain_/rotation_/stress_ 前缀区分
    """
    dir_py = Path(dir_py)
    dir_c = Path(dir_c)
    files_py = {p.name: p for p in dir_py.glob(pattern)}
    files_c = {p.name: p for p in dir_c.glob(pattern)}

    print(f"compare SAC: {dir_py}  vs  {dir_c}  ({pattern})")

    if not files_py and not files_c:
        raise AssertionError(f"no SAC files matching {pattern!r} under {dir_py} or {dir_c}")

    missing_py = sorted(set(files_c) - set(files_py))
    missing_c = sorted(set(files_py) - set(files_c))
    if missing_py:
        raise AssertionError(f"missing file(s) in {dir_py}: {missing_py}")
    if missing_c:
        raise AssertionError(f"missing file(s) in {dir_c}: {missing_c}")

    error = 0.0
    nerr = 0
    for name in sorted(files_c):
        tr_c = read(str(files_c[name]))[0]
        tr_py = read(str(files_py[name]))[0]
        if len(tr_c.data) != len(tr_py.data):
            raise AssertionError(f"npts mismatch for {name}: py={len(tr_py.data)} c={len(tr_c.data)}")
        rerr = _rel_error(tr_py.data, tr_c.data)
        if rerr is None:
            continue
        print(f"  {name}: {rerr:.6e}")
        error += rerr
        nerr += 1

    if nerr == 0:
        return 0.0
    return error / nerr


def _nc_variable_map(path: PathLike) -> dict:
    """读取 NetCDF 变量数据，跳过坐标轴 north/east/depsrc/deprcv/depth"""
    result = {}
    with netcdf_file(str(path), mmap=False) as dataset:
        for name, variable in dataset.variables.items():
            if name in {"north", "east", "depsrc", "deprcv", "depth"}:
                continue
            result[name] = np.array(variable[:], copy=True)
    return result


def compare_nc_files(path_py: PathLike, path_c: PathLike) -> float:
    """
    比较两个静态 NetCDF 文件中的物理量变量，返回平均相对误差
    """
    path_py = Path(path_py)
    path_c = Path(path_c)
    print(f"compare NC: {path_py}  vs  {path_c}")

    py_vars = _nc_variable_map(path_py)
    c_vars = _nc_variable_map(path_c)

    keys = sorted(set(py_vars) | set(c_vars))
    missing_py = sorted(set(c_vars) - set(py_vars))
    missing_c = sorted(set(py_vars) - set(c_vars))
    if missing_py:
        raise AssertionError(f"missing in py NC: {missing_py}")
    if missing_c:
        raise AssertionError(f"missing in c NC: {missing_c}")

    error = 0.0
    nerr = 0
    for key in keys:
        val_py = py_vars[key]
        val_c = c_vars[key]
        if val_py.shape != val_c.shape:
            raise AssertionError(f"shape mismatch for {key}: py={val_py.shape} c={val_c.shape}")
        rerr = _rel_error(val_py, val_c)
        if rerr is None:
            continue
        print(f"  {key}: {rerr:.6e}")
        error += rerr
        nerr += 1

    if nerr == 0:
        return 0.0
    return error / nerr


def assert_command_has(command: Sequence[object], *tokens: str) -> None:
    """断言命令列表包含给定 token"""
    values = [str(item) for item in command]
    for token in tokens:
        if token not in values:
            raise AssertionError(f"expected token {token!r} in command:\n  {' '.join(values)}")


def assert_command_equals(command: Sequence[object], expected: Iterable[str]) -> None:
    """断言命令列表与期望完全一致"""
    actual = [str(item) for item in command]
    expect = [str(item) for item in expected]
    if actual != expect:
        raise AssertionError("command mismatch:\n" f"  actual:   {' '.join(actual)}\n" f"  expected: {' '.join(expect)}")


def summarize_errors(name: str, errors: Sequence[float], threshold: float) -> float:
    """打印误差统计，超阈值则抛错"""
    arr = np.asarray(errors, dtype=float)
    mean = float(np.mean(arr)) if arr.size else 0.0
    print(f"---------------- {name} --------------------")
    print(arr)
    if arr.size:
        print(f"mean={mean:.6e}  min={np.min(arr):.6e}  max={np.max(arr):.6e}")
    if mean > threshold:
        raise AssertionError(f"{name} mean relative error {mean:.6e} exceeds {threshold:.6e}")
    return mean
