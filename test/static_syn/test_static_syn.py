import shutil
from pathlib import Path

import pygrt
from scipy.io import netcdf_file

depsrc = 2.0
deprcv = 0.0
dists = [0.0, 1.0, 2.0, 4.0, 8.0, 12.0]
norths = [-2.0, 2.0, 0.5]
easts = [-1.0, 1.0, 0.5]
modname = "../milrow"

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
    output_path="stsyn_ff.nc",
    finite_fault=ff,
    subfault_size=(1.0, 1.0),
    norths=norths,
    easts=easts,
)

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
