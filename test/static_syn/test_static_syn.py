import shutil
from pathlib import Path

import numpy as np
import pygrt
from scipy.io import netcdf_file

depsrc = 2.0
deprcv = 0.0
dists = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0]
norths = [-2.0, 2.0, 0.5]
easts = [-1.0, 1.0, 0.5]
modname = "../milrow"


def read_static_fields(path):
    with netcdf_file(path, mmap=False) as dataset:
        variables = {
            name: np.array(variable.data, dtype=np.float64, copy=True)
            for name, variable in dataset.variables.items()
        }
        attributes = {
            name: getattr(dataset, name)
            for name in ("rot2ZNE", "calc_upar")
            if hasattr(dataset, name)
        }
    return variables, attributes


def assert_zne_zrt_equivalent(zne_path, zrt_path):
    zne, zne_attributes = read_static_fields(zne_path)
    zrt, zrt_attributes = read_static_fields(zrt_path)
    assert zne_attributes["rot2ZNE"] == 1
    assert zrt_attributes["rot2ZNE"] == 0
    assert zne_attributes["calc_upar"] == 1
    assert zrt_attributes["calc_upar"] == 1
    np.testing.assert_allclose(zne["north"], zrt["north"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(zne["east"], zrt["east"], rtol=0.0, atol=0.0)

    north = zrt["north"][:, None]
    east = zrt["east"][None, :]
    distance = np.hypot(north, east)
    theta = np.where(distance == 0.0, 0.0, np.arctan2(east, north))
    cosine = np.cos(theta)
    sine = np.sin(theta)

    np.testing.assert_allclose(zne["Z"], zrt["Z"], rtol=1.0e-12, atol=1.0e-12)
    np.testing.assert_allclose(
        zne["N"], zrt["R"] * cosine - zrt["T"] * sine,
        rtol=1.0e-12, atol=1.0e-12,
    )
    np.testing.assert_allclose(
        zne["E"], zrt["R"] * sine + zrt["T"] * cosine,
        rtol=1.0e-12, atol=1.0e-12,
    )

    r = distance * 1.0e5
    s00, s01, s02 = zrt["zZ"], zrt["zR"], zrt["zT"]
    s10, s11, s12 = zrt["rZ"], zrt["rR"], zrt["rT"]
    s20, s21, s22 = zrt["tZ"], zrt["tR"], zrt["tT"]
    ur_over_r = np.array(s11, copy=True)
    ut_over_r = np.array(s12, copy=True)
    nonzero = r != 0.0
    ur_over_r[nonzero] = zrt["R"][nonzero] / r[nonzero]
    ut_over_r[nonzero] = zrt["T"][nonzero] / r[nonzero]
    converted = {
        "zZ": s00,
        "zN": s01 * cosine - s02 * sine,
        "zE": s01 * sine + s02 * cosine,
        "nZ": s10 * cosine - s20 * sine,
        "eZ": s10 * sine + s20 * cosine,
        "nN": s11 * cosine**2 + s22 * sine**2 - (s12 + s21) * sine * cosine
        + ur_over_r * sine**2 + ut_over_r * sine * cosine,
        "nE": s12 * cosine**2 - s21 * sine**2 + (s11 - s22) * sine * cosine
        - ur_over_r * sine * cosine + ut_over_r * sine**2,
        "eN": s21 * cosine**2 - s12 * sine**2 + (s11 - s22) * sine * cosine
        - ur_over_r * sine * cosine - ut_over_r * cosine**2,
        "eE": s22 * cosine**2 + s11 * sine**2 + (s12 + s21) * sine * cosine
        + ur_over_r * cosine**2 - ut_over_r * sine * cosine,
    }
    for name, expected in converted.items():
        np.testing.assert_allclose(zne[name], expected, rtol=2.0e-11, atol=1.0e-12)

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=modname)
pymod.compute_static_grn(depsrc=depsrc, deprcv=deprcv, dists=dists, calc_upar=True)

pymod.compute_static_syn(scale=1e20, output_path="stsyn.nc", norths=norths, easts=easts)
pymod.compute_static_syn(scale=1e16, output_path="stsyn.nc", force=(-1, 2, -4), norths=norths, easts=easts)
pymod.compute_static_syn(scale=1e20, output_path="stsyn.nc", strike=33, dip=44, rake=55, norths=norths, easts=easts)
pymod.compute_static_syn(scale=1e20, output_path="stsyn.nc", strike=33, dip=44, norths=norths, easts=easts)
pymod.compute_static_syn(scale=1e20, output_path="stsyn.nc", moment_tensor=(1, -2, -5, 0.5, 3, 1.2), norths=norths, easts=easts)

# ZNE / 空间导数 / 二维接收网格
pymod.compute_static_syn(scale=1e20, output_path="stsyn.nc", norths=norths, easts=easts, zne=True, calc_upar=True)
pymod.compute_static_syn(scale=1e20, output_path="stsyn_single_explicit.nc", depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)

try:
    pymod.compute_static_syn(scale=1e20, output_path="stsyn_bad.nc", depsrc=depsrc + 0.1)
    raise AssertionError("wrong single source depth should fail")
except RuntimeError:
    pass

try:
    pymod.compute_static_syn(scale=1e20, output_path="stsyn_bad.nc", deprcv=deprcv + 0.1)
    raise AssertionError("wrong single receiver depth should fail")
except RuntimeError:
    pass

# -------------------- -R 建库后合成 --------------------
pymod_r = pygrt.PyModel1D(stgrn="stgrn_r.nc", modelpath=modname)
pymod_r.compute_static_grn(depsrc=depsrc, deprcv=deprcv, dists=[0.0, 1.0, 2.0, 4.0, 8.0], calc_upar=True)
pymod_r.compute_static_syn(scale=1e20, output_path="stsyn_r.nc", norths=norths, easts=easts)

# -------------------- 多深度：点源 -Ds、深度插值、-Q、有限断层 --------------------
pymod_m = pygrt.PyModel1D(stgrn="stgrn_md.nc", modelpath=modname)
pymod_m.compute_static_grn(depsrc=[1.0, 2.0, 3.0], deprcv=0.0, dists=dists, calc_upar=True)

# Compare the two CLI outputs generated above before exercising the Python API
assert_zne_zrt_equivalent("stsyn_ff_zne_cli.nc", "stsyn_ff_zrt_cli.nc")

pymod_m.compute_static_syn(scale=1e16, output_path="stsyn_md.nc", depsrc=2.0, norths=norths, easts=easts, scale_with_mu=True)
pymod_m.compute_static_syn(
    scale=1e16, output_path="stsyn_interp.nc", depsrc=1.5,
    norths=[-2.0, 2.0, 1.0], easts=[-2.0, 2.0, 1.0],
    scale_with_mu=True, calc_upar=True,
)

rcv = Path("rcv_pts.txt")
rcv.write_text("# north east depth (km)\n0 0 0\n1 2 0\n-1 1 0\n")
pymod_m.compute_static_syn(scale=1e16, output_path="stsyn_q.nc", depsrc=2.0, recv_points=rcv)
with netcdf_file("stsyn_q.nc", mmap=False) as f:
    assert "point" in f.dimensions
    assert f.dimensions["point"] == 3

ff = Path("cfaults_tiny.inp")
# W=(2.8-1.2)/sin(90°)=1.6 km，dW=1 → 末块短于 dW，覆盖余数子断层中心
pymod_m.compute_static_syn(
    output_path="stsyn_ff_zrt.nc",
    finite_fault=ff,
    subfault_size=(1.0, 1.0),
    norths=norths,
    easts=easts,
    calc_upar=True,
)
pymod_m.compute_static_syn(
    output_path="stsyn_ff_zne.nc",
    finite_fault=ff,
    subfault_size=(1.0, 1.0),
    norths=norths,
    easts=easts,
    zne=True,
    calc_upar=True,
)
assert_zne_zrt_equivalent("stsyn_ff_zne.nc", "stsyn_ff_zrt.nc")

# The shared Coulomb fixtures cover every supported Kode and the .inr rake
# format.  Besides checking that the files are accepted, require a finite
# nonzero result so that a silently skipped source cannot pass the test.
coulomb_kodes = Path("cfaults_kodes.inp")
coulomb_rake = Path("cfaults_rake.inr")
for output_path, finite_fault in [
    ("stsyn_ff_kodes.nc", coulomb_kodes),
    ("stsyn_ff_rake.nc", coulomb_rake),
]:
    pymod_m.compute_static_syn(
        output_path=output_path,
        finite_fault=finite_fault,
        subfault_size=(1.0, 1.0),
        norths=(-4.0, 4.0, 2.0),
        easts=(-4.0, 4.0, 2.0),
    )
    with netcdf_file(output_path, mmap=False) as f:
        z = f.variables["Z"].data
        assert (z == z).all()
        assert abs(z).max() > 0.0

# 多台站深度必须给 deprcv
pymod_mr = pygrt.PyModel1D(stgrn="stgrn_mr.nc", modelpath=modname)
pymod_mr.compute_static_grn(depsrc=2.0, deprcv=[0.0, 0.5], dists=[0.0, 5.0], calc_upar=True)
pymod_mr.compute_static_syn(scale=1e20, output_path="stsyn_dr.nc", deprcv=0.25, norths=[-2.0, 2.0, 1.0], easts=[-2.0, 2.0, 1.0])

# -------------------- 错误参数 --------------------
try:
    pymod_m.compute_static_syn(scale=1e20, output_path="stsyn_bad.nc", recv_points=rcv, norths=[-1.0, 1.0, 1.0], easts=[-1.0, 1.0, 1.0])
    raise AssertionError("recv_points with norths/easts should raise")
except ValueError:
    pass

try:
    pymod_m.compute_static_syn(scale=1e20, output_path="stsyn_bad.nc", recv_points=rcv, deprcv=0.0)
    raise AssertionError("recv_points with deprcv should raise")
except ValueError:
    pass

try:
    pymod_m.compute_static_syn(scale=1e20, output_path="stsyn_bad.nc", finite_fault=ff)
    raise AssertionError("finite_fault with scale should raise")
except ValueError:
    pass

bad_dip = Path("cfaults_bad_dip.inp")
try:
    pymod_m.compute_static_syn(output_path="stsyn_bad.nc", finite_fault=bad_dip)
    raise AssertionError("dip=0 should fail in C")
except RuntimeError:
    pass

bad_bot = Path("cfaults_bad_bot.inp")
try:
    pymod_m.compute_static_syn(output_path="stsyn_bad.nc", finite_fault=bad_bot)
    raise AssertionError("bot < top should fail in C")
except RuntimeError:
    pass

for name in [
    "stgrn.nc", "stgrn_r.nc", "stgrn_md.nc", "stgrn_mr.nc",
    "stsyn.nc", "stsyn_single_explicit.nc", "stsyn_r.nc",
    "stsyn_md.nc", "stsyn_interp.nc", "stsyn_q.nc", "stsyn_ff.nc",
    "stsyn_ff_zrt.nc", "stsyn_ff_zne.nc", "stsyn_ff_zrt_cli.nc", "stsyn_ff_zne_cli.nc",
    "stsyn_ff_kodes.nc", "stsyn_ff_rake.nc", "stsyn_dr.nc", "ff_kodes_q.nc",
    "rcv_pts.txt", "cfaults_tiny.inp", "cfaults_kodes.inp", "cfaults_rake.inr",
    "cfaults_bad_dip.inp", "cfaults_bad_bot.inp",
]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)

print("test_static_syn.py: all checks passed")
