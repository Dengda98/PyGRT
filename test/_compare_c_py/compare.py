"""
端到端比较：同一组物理参数下，直接调用 grt CLI 与 Python API 的结果是否一致

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

import pygrt
from pygrt.cli import find_grt, format_float, run_grt
from compare_func import compare_nc_files, compare_sac_dirs, summarize_errors


MODEL = (HERE.parent / "milrow").resolve()
WORKDIR = HERE / "_work_compare"

DIST = 10.0
DEPSRC = 2.0
DEPRCV = 3.3
NT = 1024
DT = 0.01
SCALE = 1e24
AZ = 39.2
FN, FE, FZ = 2.0, -1.0, 4.0
STK, DIP, RAK = 77.0, 88.0, 99.0
MT = (1.0, -2.0, -5.0, 0.5, 3.0, 1.2)
NORTHS = (-3.1, 3.1, 0.6)
EASTS = (-4.1, 4.1, 0.8)

# Python 包装会显式写出默认选项；C 直调侧使用相同命令，保证参数解析一致
DYNAMIC_THRESH = 1e-10
STATIC_THRESH = 1e-10


def _clean_workdir() -> None:
    if WORKDIR.exists():
        shutil.rmtree(WORKDIR)
    WORKDIR.mkdir(parents=True)


def _grn_subdir(root: Path) -> Path:
    name = (
        f"{MODEL.name}_{format_float(DEPSRC)}_"
        f"{format_float(DEPRCV)}_{format_float(DIST)}"
    )
    return root / name


def run_c_dynamic(c_root: Path) -> None:
    """用与 Python 默认映射一致的命令行计算动态结果"""
    grn_root = c_root / "GRN"
    grn_root.mkdir(parents=True)

    run_grt([
        "greenfn",
        f"-M{MODEL}",
        f"-D{format_float(DEPSRC)}/{format_float(DEPRCV)}",
        f"-N{NT}/{format_float(DT)}+w0.8+n1",
        f"-R{format_float(DIST)}",
        f"-O{grn_root}",
        "-BfH",
        "-H-1/-1",
        "-L0",
        "-K+k50+s2+e-1",
        "-E0/0",
        "-e",
    ])

    gdir = _grn_subdir(grn_root)
    cases = [
        ("syn_ex", [], "t/0.2/0.2/0.4"),
        ("syn_sf", [f"-F{format_float(FN)}/{format_float(FE)}/{format_float(FZ)}"], "t/0.1/0.3/0.6"),
        ("syn_dc", [f"-M{format_float(STK)}/{format_float(DIP)}/{format_float(RAK)}"], "p/0.6"),
        ("syn_ts", [f"-M{format_float(STK)}/{format_float(DIP)}"], "p/0.6"),
        (
            "syn_mt",
            ["-T" + "/".join(format_float(v) for v in MT)],
            "r/3",
        ),
    ]
    for zne in (False, True):
        suffix = "-N" if zne else ""
        for name, extra, tf in cases:
            out = c_root / f"{name}{suffix}"
            cmd = [
                "syn",
                f"-G{gdir}",
                f"-A{format_float(AZ)}",
                f"-S{format_float(SCALE)}",
                f"-O{out}",
                *extra,
                f"-D{tf}",
                "-e",
            ]
            if zne:
                cmd.append("-N")
            run_grt(cmd)
            run_grt(["strain", str(out)])
            run_grt(["rotation", str(out)])
            run_grt(["stress", str(out)])


def run_py_dynamic(py_root: Path) -> None:
    """Python API 计算动态结果"""
    model = pygrt.PyModel1D(grn=py_root / "GRN", modelpath=MODEL, topbound="free", botbound="halfspace")
    model.compute_grn(
        depsrc=DEPSRC,
        deprcv=DEPRCV,
        dists=DIST,
        nt=NT,
        dt=DT,
        calc_upar=True,
        print_log=False,
    )

    cases = [
        ("syn_ex", "EX", {}, "t/0.2/0.2/0.4"),
        ("syn_sf", "SF", {"force": (FN, FE, FZ)}, "t/0.1/0.3/0.6"),
        ("syn_dc", "DC", {"strike": STK, "dip": DIP, "rake": RAK}, "p/0.6"),
        ("syn_ts", "TS", {"strike": STK, "dip": DIP}, "p/0.6"),
        ("syn_mt", "MT", {"moment_tensor": MT}, "r/3"),
    ]
    for zne in (False, True):
        suffix = "-N" if zne else ""
        for name, source, kwargs, tf in cases:
            out = py_root / f"{name}{suffix}"
            model.compute_syn(
                dist=DIST,
                azimuth=AZ,
                scale=SCALE,
                output_path=out,
                source=source,
                time_function=tf,
                zne=zne,
                calc_upar=True,
                **kwargs,
            )
            pygrt.utils.compute_strain(out)
            pygrt.utils.compute_rotation(out)
            pygrt.utils.compute_stress(out)


def run_c_static(c_root: Path) -> None:
    """用与 Python 默认映射一致的命令行计算静态结果"""
    static_dir = c_root / "static"
    static_dir.mkdir(parents=True)
    grn = static_dir / "stgrn.nc"

    run_grt([
        "static",
        "greenfn",
        f"-M{MODEL}",
        f"-D{format_float(DEPSRC)}/{format_float(DEPRCV)}",
        f"-O{grn}",
        "-BfH",
        f"-X{NORTHS[0]}/{NORTHS[1]}/{NORTHS[2]}",
        f"-Y{EASTS[0]}/{EASTS[1]}/{EASTS[2]}",
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
                *extra,
                "-e",
            ]
            if zne:
                cmd.append("-N")
            run_grt(cmd)
            run_grt(["static", "strain", str(out)])
            run_grt(["static", "rotation", str(out)])
            run_grt(["static", "stress", str(out)])


def run_py_static(py_root: Path) -> None:
    """Python API 计算静态结果"""
    static_dir = py_root / "static"
    static_dir.mkdir(parents=True)
    model = pygrt.PyModel1D(stgrn=static_dir / "stgrn.nc", modelpath=MODEL)
    model.compute_static_grn(
        depsrc=DEPSRC,
        deprcv=DEPRCV,
        norths=NORTHS,
        easts=EASTS,
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
                zne=zne,
                calc_upar=True,
                **kwargs,
            )
            pygrt.utils.compute_strain(out)
            pygrt.utils.compute_rotation(out)
            pygrt.utils.compute_stress(out)


def compare_dynamic(c_root: Path, py_root: Path) -> list[float]:
    errors = []
    # 格林函数
    errors.append(
        compare_sac_dirs(_grn_subdir(py_root / "GRN"), _grn_subdir(c_root / "GRN"))
    )

    for zne in (False, True):
        suffix = "-N" if zne else ""
        for name in ("syn_ex", "syn_sf", "syn_dc", "syn_ts", "syn_mt"):
            py_dir = py_root / f"{name}{suffix}"
            c_dir = c_root / f"{name}{suffix}"
            # 位移三分量
            disp = "[ZNE].sac" if zne else "[ZRT].sac"
            errors.append(compare_sac_dirs(py_dir, c_dir, pattern=disp))
            # 空间导数（两位：小写方向 + 大写分量）
            deriv = "[zne][ZNE].sac" if zne else "[zrt][ZRT].sac"
            errors.append(compare_sac_dirs(py_dir, c_dir, pattern=deriv))
            # 张量按文件名前缀分别比对，避免 channel 撞名漏比
            for prefix in ("strain_", "rotation_", "stress_"):
                errors.append(compare_sac_dirs(py_dir, c_dir, pattern=f"{prefix}*.sac"))
    return errors


def compare_static(c_root: Path, py_root: Path) -> list[float]:
    errors = []
    errors.append(
        compare_nc_files(py_root / "static" / "stgrn.nc", c_root / "static" / "stgrn.nc")
    )
    for zne in (False, True):
        suffix = "-N" if zne else ""
        for name in ("stsyn_ex", "stsyn_sf", "stsyn_dc", "stsyn_ts", "stsyn_mt"):
            errors.append(
                compare_nc_files(
                    py_root / "static" / f"{name}{suffix}.nc",
                    c_root / "static" / f"{name}{suffix}.nc",
                )
            )
    return errors


def main():
    print(f"using grt: {find_grt()}")
    _clean_workdir()
    c_root = WORKDIR / "c"
    py_root = WORKDIR / "py"
    c_root.mkdir()
    py_root.mkdir()

    print("=== dynamic: C CLI ===")
    run_c_dynamic(c_root)
    print("=== dynamic: Python API ===")
    run_py_dynamic(py_root)
    dyn_errors = compare_dynamic(c_root, py_root)

    print("=== static: C CLI ===")
    run_c_static(c_root)
    print("=== static: Python API ===")
    run_py_static(py_root)
    st_errors = compare_static(c_root, py_root)

    summarize_errors("dynamic", dyn_errors, DYNAMIC_THRESH)
    summarize_errors("static", st_errors, STATIC_THRESH)
    print("All end-to-end comparisons passed.")


if __name__ == "__main__":
    main()
