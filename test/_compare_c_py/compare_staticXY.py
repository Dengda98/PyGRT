"""
静态合成时指定新 XY 网格：比较 C CLI 与 Python API 结果

两侧最终都执行同一 grt 可执行文件，因此本测试主要验证 Python 参数拼装与
手写 C 命令是否一致，而非两套独立数值实现之间的交叉校验
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import numpy as np
import pygrt
from pygrt.cli import find_grt, format_float, run_grt
from compare_func import compare_nc_files, summarize_errors


MODEL = (HERE.parent / "milrow").resolve()
WORKDIR = HERE / "_work_compare_staticXY"

DEPSRC = 2.0
DEPRCV = 3.3
SCALE = 1e24
FN, FE, FZ = 2.0, -1.0, 4.0
STK, DIP, RAK = 77.0, 88.0, 99.0
MT = (1.0, -2.0, -5.0, 0.5, 3.0, 1.2)
# 格林函数库用震中距采样；合成时再插值到新网格
DISTARR = np.arange(0.0, 10.0 + 1e-8, 0.1)
NORTHS = (-3.1, 3.1, 0.6)
EASTS = (-4.1, 4.1, 0.8)
THRESH = 1e-10


def _clean_workdir() -> None:
    if WORKDIR.exists():
        shutil.rmtree(WORKDIR)
    WORKDIR.mkdir(parents=True)


def run_c(c_root: Path) -> None:
    static_dir = c_root / "static"
    static_dir.mkdir(parents=True)
    grn = static_dir / "stgrn.nc"
    rs = ",".join(format_float(float(r)) for r in DISTARR)

    run_grt([
        "static",
        "greenfn",
        f"-M{MODEL}",
        f"-D{format_float(DEPSRC)}/{format_float(DEPRCV)}",
        f"-O{grn}",
        "-BfH",
        f"-R{rs}",
        "-L15",
        "-K+k50+e-1",
        "-e",
    ])

    cases = [
        ("stsyn_ex", []),
        ("stsyn_sf", [f"-F{format_float(FN)}/{format_float(FE)}/{format_float(FZ)}"]),
        ("stsyn_dc", [f"-M{format_float(STK)}/{format_float(DIP)}/{format_float(RAK)}"]),
        ("stsyn_ts", [f"-M{format_float(STK)}/{format_float(DIP)}"]),
        ("stsyn_mt", ["-T" + "/".join(format_float(v) for v in MT)]),
    ]
    for zne in (False, True):
        suffix = "-N" if zne else ""
        for name, extra in cases:
            out = static_dir / f"{name}{suffix}.nc"
            cmd = [
                "static",
                "syn",
                f"-G{grn}",
                f"-S{format_float(SCALE)}",
                f"-O{out}",
                f"-X{NORTHS[0]}/{NORTHS[1]}/{NORTHS[2]}",
                f"-Y{EASTS[0]}/{EASTS[1]}/{EASTS[2]}",
                *extra,
                "-e",
            ]
            if zne:
                cmd.append("-N")
            run_grt(cmd)
            run_grt(["static", "strain", str(out)])
            run_grt(["static", "rotation", str(out)])
            run_grt(["static", "stress", str(out)])


def run_py(py_root: Path) -> None:
    static_dir = py_root / "static"
    static_dir.mkdir(parents=True)
    model = pygrt.PyModel1D(stgrn=static_dir / "stgrn.nc", modelpath=MODEL)
    model.compute_static_grn(
        depsrc=DEPSRC,
        deprcv=DEPRCV,
        distarr=DISTARR,
        calc_upar=True,
    )

    cases = [
        ("stsyn_ex", "EX", {}),
        ("stsyn_sf", "SF", {"force": (FN, FE, FZ)}),
        ("stsyn_dc", "DC", {"strike": STK, "dip": DIP, "rake": RAK}),
        ("stsyn_ts", "TS", {"strike": STK, "dip": DIP}),
        ("stsyn_mt", "MT", {"moment_tensor": MT}),
    ]
    for zne in (False, True):
        suffix = "-N" if zne else ""
        for name, source, kwargs in cases:
            out = static_dir / f"{name}{suffix}.nc"
            model.compute_static_syn(
                scale=SCALE,
                output_path=out,
                source=source,
                norths=NORTHS,
                easts=EASTS,
                zne=zne,
                calc_upar=True,
                **kwargs,
            )
            pygrt.utils.compute_strain(out)
            pygrt.utils.compute_rotation(out)
            pygrt.utils.compute_stress(out)


def main():
    print(f"using grt: {find_grt()}")
    _clean_workdir()
    c_root = WORKDIR / "c"
    py_root = WORKDIR / "py"
    c_root.mkdir()
    py_root.mkdir()

    print("=== staticXY: C CLI ===")
    run_c(c_root)
    print("=== staticXY: Python API ===")
    run_py(py_root)

    errors = [
        compare_nc_files(
            py_root / "static" / "stgrn.nc",
            c_root / "static" / "stgrn.nc",
        )
    ]
    for zne in (False, True):
        suffix = "-N" if zne else ""
        for name in ("stsyn_ex", "stsyn_sf", "stsyn_dc", "stsyn_ts", "stsyn_mt"):
            errors.append(
                compare_nc_files(
                    py_root / "static" / f"{name}{suffix}.nc",
                    c_root / "static" / f"{name}{suffix}.nc",
                )
            )

    summarize_errors("staticXY", errors, THRESH)
    print("All staticXY comparisons passed.")


if __name__ == "__main__":
    main()
