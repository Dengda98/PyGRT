"""
Compare STGRNLIB nc results:
1) CLI vs Python for single/multi source and receiver depths
2) Multi-depth library slices vs corresponding single-depth results
3) Synthesis from multi-depth libraries: CLI vs Python, exact-depth vs
   single-depth library, depth-interpolation linearity, and tensor postprocess
"""

from pathlib import Path

import numpy as np
import pygrt
from pygrt.cli import format_float, run_grt
from scipy.io import netcdf_file

from compare_func import compare_nc_files

# 与 shell 脚本保持一致
NORTHS = (-2.0, 2.0, 1.0)
EASTS = (-2.0, 2.0, 1.0)
DEPSRCS = np.array([1.0, 2.0, 3.0])
DEPRCVS = np.array([0.0, 0.5, 1.0])

MODNAME = str((Path(__file__).resolve().parent.parent / "milrow"))
CMPDIR = Path("stgrnlib_cmp")
ATOL = 1e-10
RTOL = 1e-8
# 平均相对误差阈值
MEAN_RERR_MAX = 1e-6
SCALE = 1e20

SYN_SKIP = {
    "north", "east", "depth", "depsrc", "deprcv",
    "rcv_va", "rcv_vb", "rcv_rho", "model",
}


def _coord_tag(z: float) -> str:
    """与 shell 脚本一致：整数写 0/1，小数点换成 p（如 0p5）"""
    zf = float(z)
    if zf.is_integer():
        return str(int(zf))
    return str(zf).replace(".", "p")


def load_stgrnlib_nc(path: str) -> dict:
    """读取 STGRNLIB nc 为字典（坐标 + 各通道数组）"""
    out = {}
    with netcdf_file(path, mmap=False) as f:
        for k, v in f.variables.items():
            out[k] = np.array(v[:], dtype=float).copy()
    return out


def compare_nc_data(a: dict, b: dict, label: str) -> float:
    """逐通道比较，返回平均相对误差"""
    for key in ("depsrc", "deprcv", "north", "east"):
        if not np.allclose(a[key], b[key], rtol=RTOL, atol=ATOL):
            raise ValueError(f"{label}: coordinate mismatch on '{key}'\n"
                             f"  a={a[key]}\n  b={b[key]}")

    med_keys = ("src_va", "src_vb", "src_rho", "rcv_va", "rcv_vb", "rcv_rho")
    for key in med_keys:
        if key in a and key in b:
            if not np.allclose(a[key], b[key], rtol=RTOL, atol=ATOL):
                raise ValueError(f"{label}: medium mismatch on '{key}'\n"
                                 f"  a={a[key]}\n  b={b[key]}")

    skip = {"depsrc", "deprcv", "north", "east", "model", *med_keys}
    keys = sorted(set(a.keys()) & set(b.keys()) - skip)
    if not keys:
        raise ValueError(f"{label}: no channel variables to compare")

    errors = []
    for k in keys:
        va, vb = a[k], b[k]
        if va.ndim != 4 or vb.ndim != 4:
            continue
        if va.shape != vb.shape:
            raise ValueError(f"{label}: shape mismatch on '{k}': {va.shape} vs {vb.shape}")
        if np.all(va == 0.0) and np.all(vb == 0.0):
            continue
        denom = np.mean(np.abs(vb))
        if denom == 0.0:
            denom = np.mean(np.abs(va))
        rerr = np.sum(np.abs(va - vb)) / denom
        print(f"  {k}: {rerr:.6e}")
        errors.append(rerr)

    mean_err = float(np.mean(errors)) if errors else 0.0
    print(f"{label}: mean relative error = {mean_err:.6e}")
    if mean_err > MEAN_RERR_MAX:
        raise ValueError(f"{label}: mean relative error too large ({mean_err})")
    return mean_err


def extract_slice(lib: dict, isrc: int, ircv: int) -> dict:
    """从多深度库取出单一 (depsrc, deprcv) 切片，维度变为 [1,1,north,east]"""
    out = {
        "depsrc": np.array([lib["depsrc"][isrc]], dtype=float),
        "deprcv": np.array([lib["deprcv"][ircv]], dtype=float),
        "north": lib["north"].copy(),
        "east": lib["east"].copy(),
        "src_va": np.array([lib["src_va"][isrc]], dtype=float),
        "src_vb": np.array([lib["src_vb"][isrc]], dtype=float),
        "src_rho": np.array([lib["src_rho"][isrc]], dtype=float),
        "rcv_va": np.array([lib["rcv_va"][ircv]], dtype=float),
        "rcv_vb": np.array([lib["rcv_vb"][ircv]], dtype=float),
        "rcv_rho": np.array([lib["rcv_rho"][ircv]], dtype=float),
    }
    skip = set(out.keys()) | {"model"}
    for k, v in lib.items():
        if k in skip:
            continue
        if getattr(v, "ndim", 0) != 4:
            continue
        # 通道数据 shape: (ndepsrc, ndeprcv, nnorth, neast)
        out[k] = v[isrc:isrc + 1, ircv:ircv + 1, :, :].copy()
    return out


def py_compute_to_nc(depsrcs, deprcvs, outpath: str):
    pymod = pygrt.PyModel1D(MODNAME)
    pymod.set_static_grn_path(outpath)
    pymod.compute_static_grn(
        depsrc=depsrcs, deprcv=deprcvs, norths=NORTHS, easts=EASTS, calc_upar=True,
    )
    assert Path(outpath).is_file()


def load_syn_fields(path: Path) -> dict:
    """读取合成 nc 中的场量（跳过坐标与介质）"""
    out = {}
    with netcdf_file(str(path), mmap=False) as f:
        for k, v in f.variables.items():
            if k in SYN_SKIP:
                continue
            out[k] = np.array(v[:], dtype=float).copy()
    return out


def compare_syn_fields(a: dict, b: dict, label: str) -> float:
    """比较两个合成场量字典，返回平均相对误差"""
    keys = sorted(set(a.keys()) & set(b.keys()))
    if not keys:
        raise ValueError(f"{label}: no syn fields to compare")
    errors = []
    for k in keys:
        va, vb = a[k], b[k]
        if va.shape != vb.shape:
            raise ValueError(f"{label}: shape mismatch on '{k}': {va.shape} vs {vb.shape}")
        if np.all(va == 0.0) and np.all(vb == 0.0):
            continue
        denom = np.mean(np.abs(vb))
        if denom == 0.0:
            denom = np.mean(np.abs(va))
        rerr = np.sum(np.abs(va - vb)) / denom
        print(f"  {k}: {rerr:.6e}")
        errors.append(rerr)
    mean_err = float(np.mean(errors)) if errors else 0.0
    print(f"{label}: mean relative error = {mean_err:.6e}")
    if mean_err > MEAN_RERR_MAX:
        raise ValueError(f"{label}: mean relative error too large ({mean_err})")
    return mean_err


def c_static_syn(grn: Path, out: Path, extra: list) -> None:
    run_grt(["static", "syn", f"-G{grn}", f"-S{format_float(SCALE)}", f"-O{out}", *extra, "-e"])


def py_static_syn(grn: Path, out: Path, **kwargs) -> None:
    model = pygrt.PyModel1D(MODNAME)
    model.set_static_grn_path(grn)
    model.compute_static_syn(scale=SCALE, output_path=out, calc_upar=True, **kwargs)


def compare_syn_cli_py() -> list:
    """同一多深度库上，CLI 与 Python 合成结果应一致"""
    errors = []
    grn = CMPDIR / "stgrn_mm.nc"

    print("\n--- syn exact depth Ds=2 Dr=0.5 ---")
    c_out = CMPDIR / "stsyn_exact_c.nc"
    py_out = CMPDIR / "stsyn_exact_py.nc"
    c_static_syn(grn, c_out, ["-Ds2", "-Dr0.5"])
    py_static_syn(grn, py_out, depsrc=2.0, deprcv=0.5)
    errors.append(compare_nc_files(py_out, c_out))

    print("\n--- syn interpolated Ds=1.5 Dr=0.25 + new XY ---")
    c_out = CMPDIR / "stsyn_interp_c.nc"
    py_out = CMPDIR / "stsyn_interp_py.nc"
    xy = ["-Ds1.5", "-Dr0.25", "-X-1/1/1", "-Y-1/1/1"]
    c_static_syn(grn, c_out, xy)
    py_static_syn(
        grn, py_out, depsrc=1.5, deprcv=0.25,
        norths=(-1.0, 1.0, 1.0), easts=(-1.0, 1.0, 1.0),
    )
    errors.append(compare_nc_files(py_out, c_out))

    print("\n--- syn -Q points ---")
    rcv = CMPDIR / "rcv_pts.txt"
    rcv.write_text("0 0 0\n1 1 0.5\n-1 0.5 1\n")
    c_out = CMPDIR / "stsyn_q_c.nc"
    py_out = CMPDIR / "stsyn_q_py.nc"
    c_static_syn(grn, c_out, ["-Ds2", f"-Q{rcv}"])
    py_static_syn(grn, py_out, depsrc=2.0, recv_points=rcv)
    errors.append(compare_nc_files(py_out, c_out))

    print("\n--- syn + strain/rotation/stress ---")
    c_ten = CMPDIR / "stsyn_ten_c.nc"
    py_ten = CMPDIR / "stsyn_ten_py.nc"
    c_static_syn(grn, c_ten, ["-Ds2", "-Dr0.5", "-N"])
    py_static_syn(grn, py_ten, depsrc=2.0, deprcv=0.5, zne=True)
    run_grt(["static", "strain", str(c_ten)])
    run_grt(["static", "rotation", str(c_ten)])
    run_grt(["static", "stress", str(c_ten)])
    pygrt.utils.compute_strain(py_ten)
    pygrt.utils.compute_rotation(py_ten)
    pygrt.utils.compute_stress(py_ten)
    errors.append(compare_nc_files(py_ten, c_ten))
    return errors


def compare_syn_multi_vs_single() -> list:
    """多深度库在节点深度上的合成，应与对应单深度库一致"""
    print("\n--- syn multi-depth lib at Ds=2 Dr=0.5 vs single-depth lib ---")
    mm_out = CMPDIR / "stsyn_mm_node.nc"
    ref_out = CMPDIR / "stsyn_ref_node.nc"
    c_static_syn(CMPDIR / "stgrn_mm.nc", mm_out, ["-Ds2", "-Dr0.5"])
    c_static_syn(CMPDIR / "stgrn_ref_zs2_zr0p5.nc", ref_out, [])
    return [compare_syn_fields(load_syn_fields(mm_out), load_syn_fields(ref_out),
                               "syn multi vs single [zs=2,zr=0.5]")]


def compare_syn_depth_linearity() -> list:
    """线性插值：syn(1.5) 应等于 0.5*(syn(1)+syn(2))（同一 deprcv，不用 -Su）"""
    print("\n--- syn depth interpolation linearity Ds=1.5 vs 0.5*(1+2) ---")
    grn = CMPDIR / "stgrn_mm.nc"
    out1 = CMPDIR / "stsyn_lin_z1.nc"
    out2 = CMPDIR / "stsyn_lin_z2.nc"
    outm = CMPDIR / "stsyn_lin_z15.nc"
    c_static_syn(grn, out1, ["-Ds1", "-Dr0"])
    c_static_syn(grn, out2, ["-Ds2", "-Dr0"])
    c_static_syn(grn, outm, ["-Ds1.5", "-Dr0"])
    a = load_syn_fields(out1)
    b = load_syn_fields(out2)
    mid = load_syn_fields(outm)
    pred = {k: 0.5 * (a[k] + b[k]) for k in mid if k in a and k in b}
    return [compare_syn_fields(pred, mid, "syn depth linearity [Ds=1.5]")]


def main():
    all_errs = []

    cases = [
        ("ss", 2.0, 0.0),
        ("ms", DEPSRCS, 0.0),
        ("mr", 2.0, DEPRCVS),
        ("mm", DEPSRCS, DEPRCVS),
    ]

    print("================ CLI vs Python ================")
    for tag, zs, zr in cases:
        c_path = CMPDIR / f"stgrn_{tag}.nc"
        py_path = CMPDIR / f"stgrn_{tag}_py.nc"
        print(f"\n--- case {tag} ---")
        py_compute_to_nc(zs, zr, str(py_path))
        c_data = load_stgrnlib_nc(str(c_path))
        py_data = load_stgrnlib_nc(str(py_path))
        all_errs.append(compare_nc_data(c_data, py_data, f"CLI vs Python [{tag}]"))

    print("\n================ Multi vs Single (CLI) ================")
    mm = load_stgrnlib_nc(str(CMPDIR / "stgrn_mm.nc"))
    for isrc, zs in enumerate(DEPSRCS):
        for ircv, zr in enumerate(DEPRCVS):
            ref_path = CMPDIR / f"stgrn_ref_zs{_coord_tag(zs)}_zr{_coord_tag(zr)}.nc"
            print(f"\n--- slice depsrc={zs}, deprcv={zr} ---")
            sl = extract_slice(mm, isrc, ircv)
            ref = load_stgrnlib_nc(str(ref_path))
            all_errs.append(compare_nc_data(sl, ref, f"multi vs single CLI [zs={zs},zr={zr}]"))

    print("\n================ Multi vs Single (Python) ================")
    mm_py_path = CMPDIR / "stgrn_mm_py.nc"
    mm_py = load_stgrnlib_nc(str(mm_py_path))
    for isrc, zs in enumerate(DEPSRCS):
        for ircv, zr in enumerate(DEPRCVS):
            single_path = CMPDIR / f"stgrn_py_zs{_coord_tag(zs)}_zr{_coord_tag(zr)}.nc"
            print(f"\n--- py slice depsrc={zs}, deprcv={zr} ---")
            py_compute_to_nc(float(zs), float(zr), str(single_path))
            sl = extract_slice(mm_py, isrc, ircv)
            ref = load_stgrnlib_nc(str(single_path))
            all_errs.append(compare_nc_data(sl, ref, f"multi vs single Python [zs={zs},zr={zr}]"))

    print("\n================ Syn: CLI vs Python ================")
    all_errs.extend(compare_syn_cli_py())

    print("\n================ Syn: multi vs single (exact depth) ================")
    all_errs.extend(compare_syn_multi_vs_single())

    print("\n================ Syn: depth interpolation linearity ================")
    all_errs.extend(compare_syn_depth_linearity())

    all_errs = np.array(all_errs, dtype=float)
    print("\n================ Summary ================")
    print(all_errs)
    print(f"mean={np.mean(all_errs):.6e}  min={np.min(all_errs):.6e}  max={np.max(all_errs):.6e}")
    if np.mean(all_errs) > MEAN_RERR_MAX:
        raise ValueError("overall mean relative error too large")
    print("compare_stgrnlib.py: all checks passed")


if __name__ == "__main__":
    main()
