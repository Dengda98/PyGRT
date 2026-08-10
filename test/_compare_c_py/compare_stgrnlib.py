"""
Compare STGRNLIB nc results:
1) CLI vs Python for single/multi source and receiver depths
2) Multi-depth library slices vs corresponding single-depth results
"""

from pathlib import Path

import numpy as np
import pygrt
from scipy.io import netcdf_file

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

    skip = {"depsrc", "deprcv", "north", "east", *med_keys}
    keys = sorted(set(a.keys()) & set(b.keys()) - skip)
    if not keys:
        raise ValueError(f"{label}: no channel variables to compare")

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
    skip = set(out.keys())
    for k, v in lib.items():
        if k in skip:
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

    all_errs = np.array(all_errs, dtype=float)
    print("\n================ Summary ================")
    print(all_errs)
    print(f"mean={np.mean(all_errs):.6e}  min={np.min(all_errs):.6e}  max={np.max(all_errs):.6e}")
    if np.mean(all_errs) > MEAN_RERR_MAX:
        raise ValueError("overall mean relative error too large")
    print("compare_stgrnlib.py: all checks passed")


if __name__ == "__main__":
    main()
