import shutil
from pathlib import Path

import numpy as np
import pygrt
from pygrt.cli import run_grt
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


def assert_finite_source_zne(reference_path, actual_path):
    reference, reference_attributes = read_static_fields(reference_path)
    actual, actual_attributes = read_static_fields(actual_path)
    assert reference_attributes == actual_attributes == {"rot2ZNE": 1, "calc_upar": 1}
    assert set(reference) == set(actual)
    assert {"Z", "N", "E", "zZ", "nN", "eE"} <= set(actual)
    assert not {"R", "T", "rR", "tT"} & set(actual)
    for name in reference:
        np.testing.assert_allclose(actual[name], reference[name], rtol=2e-11, atol=1e-12)


def assert_receiver_geometry(path):
    with netcdf_file(path, mmap=False) as dataset:
        for name in ("strike", "dip", "rake"):
            assert name in dataset.variables
            assert dataset.variables[name].data.shape == (3,)
        np.testing.assert_allclose(dataset.variables["strike"].data, [10.0, 40.0, 70.0])
        np.testing.assert_allclose(dataset.variables["dip"].data, [20.0, 50.0, 80.0])
        np.testing.assert_allclose(dataset.variables["rake"].data, [30.0, 60.0, 90.0])


def read_fault_receiver(path):
    with netcdf_file(path, mmap=False) as dataset:
        layout = getattr(dataset, "layout")
        if isinstance(layout, bytes):
            layout = layout.decode()
        assert layout == "faults"
        point = dataset.dimensions["point"]
        nfault = dataset.dimensions["nfault"]
        offset = np.array(dataset.variables["offset"].data, dtype=np.int64, copy=True)
        assert offset.shape == (nfault,)
        assert offset[-1] == point
        point_fields = {
            name: np.array(variable.data, dtype=np.float64, copy=True)
            for name, variable in dataset.variables.items()
            if variable.data.ndim == 1 and variable.data.shape[0] == point
        }
        faults = []
        for ifault in range(nfault):
            start = 0 if ifault == 0 else offset[ifault - 1]
            end = offset[ifault]
            coordinates = np.column_stack([
                point_fields["north"][start:end],
                point_fields["east"][start:end],
                point_fields["depth"][start:end],
            ])
            fields = {name: values[start:end] for name, values in point_fields.items()}
            faults.append({
                "coordinates": coordinates,
                "fields": fields,
            })
        geometry = {
            name: np.array(dataset.variables[name].data, dtype=np.float64, copy=True)
            for name in ("strike", "dip", "rake", "offset", "stksize", "dipsize")
        }
    return faults, geometry


def write_receiver_points(path, faults):
    coordinates = np.concatenate([fault["coordinates"] for fault in faults])
    lines = ["# north east depth\n"]
    lines.extend(f"{north:.15g} {east:.15g} {depth:.15g}\n" for north, east, depth in coordinates)
    path.write_text("".join(lines))
    return coordinates


def assert_faults_match_points(fault_path, points_path, field_names, expected_subsizes=(9, 12)):
    faults, geometry = read_fault_receiver(fault_path)
    assert [fault["coordinates"].shape[0] for fault in faults] == list(expected_subsizes)
    np.testing.assert_array_equal(geometry["offset"], np.cumsum(expected_subsizes))
    assert geometry["strike"].shape == (2,)
    assert geometry["dip"].shape == (2,)
    assert geometry["rake"].shape == (2,)
    assert geometry["stksize"].shape == (2,)
    assert geometry["dipsize"].shape == (2,)
    assert geometry["rake"][1] == -999.0

    with netcdf_file(points_path, mmap=False) as dataset:
        coordinates = np.column_stack([
            np.array(dataset.variables[name].data, dtype=np.float64, copy=True)
            for name in ("north", "east", "depth")
        ])
        assert coordinates.shape == (21, 3)
        np.testing.assert_allclose(coordinates, np.concatenate([fault["coordinates"] for fault in faults]))
        for name in field_names:
            expected = np.concatenate([fault["fields"][name] for fault in faults])
            np.testing.assert_allclose(
                expected, np.array(dataset.variables[name].data, dtype=np.float64, copy=True),
                rtol=1.0e-12, atol=1.0e-12,
            )


def compare_fault_postprocess(fault_path, points_path, prefix, component_pairs):
    faults, _ = read_fault_receiver(fault_path)
    with netcdf_file(points_path, mmap=False) as dataset:
        for first, second in component_pairs:
            name = f"{prefix}_{first}{second}"
            expected = np.concatenate([
                fault["fields"][name] for fault in faults
            ])
            np.testing.assert_allclose(
                expected, np.array(dataset.variables[name].data, dtype=np.float64, copy=True),
                rtol=1.0e-12, atol=1.0e-12,
            )

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=modname)
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists, calc_upar=True)

# 输出元数据应标记实际基本源型
for source_type, mechanism in [
    ("EX", {}),
    ("SF", {"force": (-1, 2, -4)}),
    ("DC", {"strike": 33, "dip": 44, "rake": 55}),
    ("TS", {"strike": 33, "dip": 44}),
    ("MT", {"moment_tensor": (1, -2, -5, 0.5, 3, 1.2)}),
]:
    pymod.static_syn(scale=1e16 if source_type == "SF" else 1e20, output_path="stsyn.nc", norths=norths, easts=easts, **mechanism)
    with netcdf_file("stsyn.nc", mmap=False) as dataset:
        assert dataset.computeType.decode() == source_type

# ZNE / 空间导数 / 二维接收网格
pymod.static_syn(scale=1e20, output_path="stsyn.nc", norths=norths, easts=easts, zne=True, calc_upar=True)
pymod.static_syn(scale=1e20, output_path="stsyn_single_explicit.nc", depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)

try:
    pymod.static_syn(scale=1e20, output_path="stsyn_bad.nc", depsrc=depsrc + 0.1)
    raise AssertionError("wrong single source depth should fail")
except RuntimeError:
    pass

try:
    pymod.static_syn(scale=1e20, output_path="stsyn_bad.nc", deprcv=deprcv + 0.1)
    raise AssertionError("wrong single receiver depth should fail")
except RuntimeError:
    pass

# 单节点限制应优先于剖分和几何检查，显式 +i 及点源 Kode 均不能绕过
pymod_node = pygrt.PyModel1D(stgrn="stgrn_node.nc", modelpath=modname)
pymod_node.static_syn(scale=1e20, output_path="stsyn_node_api.nc", calc_upar=True)
for kode in (100, 400, 500):
    for subdivision in ("", "+i2/2"):
        try:
            run_grt([
                "static_syn", "-Gstgrn_node.nc", f"-Ccfaults_node_{kode}.inp{subdivision}",
                "-Ostsyn_node_bad.nc",
            ], print_log=False)
        except RuntimeError as error:
            assert "Finite sources (-C) require a Green-function library with at least two nodes." in str(error)
        else:
            raise AssertionError("single-node library should reject every -C source")
        assert not Path("stsyn_node_bad.nc").exists()

try:
    pymod_node.static_syn(output_path="stsyn_node_bad.nc", src_fault="cfaults_node_400.inp", src_fault_size=(2.0, 2.0))
    raise AssertionError("Python API should reject source faults with a single-node library")
except RuntimeError:
    pass
assert not Path("stsyn_node_bad.nc").exists()

# 只有一个水平点而深度轴包含多个采样时，仍属于多节点库
for axis, source_depths, receiver_depths, distances in (
    ("src", [1.0, 3.0], 0.0, [1.0]),
    ("rcv", 2.0, [0.0, 1.0], [1.0]),
    ("dist", 2.0, 0.0, [1.0, 2.0]),
):
    pymod_two = pygrt.PyModel1D(stgrn=f"stgrn_two_{axis}.nc", modelpath=modname)
    pymod_two.static_greenfn(depsrc=source_depths, deprcv=receiver_depths, dists=distances)
    for kode in (400, 500):
        pymod_two.static_syn(
            output_path=f"stsyn_two_{axis}_{kode}.nc", src_fault=f"cfaults_node_{kode}.inp", deprcv=0.0,
        )

# -------------------- -R 建库后合成 --------------------
pymod_r = pygrt.PyModel1D(stgrn="stgrn_r.nc", modelpath=modname)
pymod_r.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=[0.0, 1.0, 2.0, 4.0, 8.0], calc_upar=True)
pymod_r.static_syn(scale=1e20, output_path="stsyn_r.nc", norths=norths, easts=easts)

# -------------------- 多深度：点源 -Ds、深度插值、-Q、有限断层 --------------------
pymod_m = pygrt.PyModel1D(stgrn="stgrn_md.nc", modelpath=modname)
pymod_m.static_greenfn(depsrc=[1.0, 2.0, 3.0], deprcv=0.0, dists=dists, calc_upar=True)

# Compare the two CLI outputs generated above before exercising the Python API
assert_finite_source_zne("stsyn_ff_zne_cli.nc", "stsyn_ff_zrt_cli.nc")
assert_receiver_geometry("stsyn_q6.nc")

syn_field_names = (
    "Z", "R", "T", "zZ", "zR", "zT", "rZ", "rR", "rT", "tZ", "tR", "tT",
)
faults, fault_geometry = read_fault_receiver("stsyn_rf.nc")
np.testing.assert_array_equal(fault_geometry["stksize"], [3, 3])
np.testing.assert_array_equal(fault_geometry["dipsize"], [3, 4])
default_faults, default_geometry = read_fault_receiver("stsyn_rf_default.nc")
assert [fault["coordinates"].shape[0] for fault in default_faults] == [1, 1]
np.testing.assert_array_equal(default_geometry["offset"], [1, 2])
np.testing.assert_array_equal(default_geometry["stksize"], [1, 1])
np.testing.assert_array_equal(default_geometry["dipsize"], [1, 1])
assert default_geometry["rake"][1] == -999.0
rcv_fault_q = Path("rcv_faults_q.txt")
write_receiver_points(rcv_fault_q, faults)
run_grt(
    ["static_syn", "-Gstgrn_rf.nc", "-Su1e16", "-Ds2", f"-Q{rcv_fault_q}", "-e", "-Ostsyn_rq.nc"],
    print_log=False,
)
assert_faults_match_points("stsyn_rf.nc", "stsyn_rq.nc", syn_field_names)

# Python API 的 -R 路径与命令行输出保持一致
pymod_rf = pygrt.PyModel1D(stgrn="stgrn_rf.nc", modelpath=modname)
pymod_rf.static_syn(
    scale=1e16, output_path="stsyn_rf_api.nc", depsrc=2.0,
    rcv_fault="rcv_faults.inp", rcv_fault_size=(0.75, 0.75), scale_with_mu=True,
    calc_upar=True,
)
api_faults, _ = read_fault_receiver("stsyn_rf_api.nc")
for shell_fault, api_fault in zip(faults, api_faults):
    np.testing.assert_allclose(shell_fault["coordinates"], api_fault["coordinates"])
    for name in syn_field_names:
        np.testing.assert_allclose(shell_fault["fields"][name], api_fault["fields"][name])

for module, component_pairs in (
    ("strain", (("Z", "Z"), ("Z", "R"), ("Z", "T"), ("R", "R"), ("R", "T"), ("T", "T"))),
    ("stress", (("Z", "Z"), ("Z", "R"), ("Z", "T"), ("R", "R"), ("R", "T"), ("T", "T"))),
    ("rotation", (("Z", "R"), ("Z", "T"), ("R", "T"))),
):
    run_grt([f"static_{module}", "stsyn_rf.nc"], print_log=False)
    run_grt([f"static_{module}", "stsyn_rq.nc"], print_log=False)
    compare_fault_postprocess("stsyn_rf.nc", "stsyn_rq.nc", module, component_pairs)

pymod_m.static_syn(scale=1e16, output_path="stsyn_md.nc", depsrc=2.0, norths=norths, easts=easts, scale_with_mu=True)
pymod_m.static_syn(
    scale=1e16, output_path="stsyn_interp.nc", depsrc=1.5,
    norths=[-2.0, 2.0, 1.0], easts=[-2.0, 2.0, 1.0],
    scale_with_mu=True, calc_upar=True,
)

rcv = Path("rcv_pts.txt")
rcv.write_text("# north east depth (km)\n0 0 0\n1 2 0\n-1 1 0\n")
pymod_m.static_syn(scale=1e16, output_path="stsyn_q.nc", depsrc=2.0, rcv_points=rcv)
with netcdf_file("stsyn_q.nc", mmap=False) as f:
    assert "point" in f.dimensions
    assert f.dimensions["point"] == 3

ff = Path("cfaults_tiny.inp")
# W=(2.8-1.2)/sin(90°)=1.6 km，dW=1 → 末块短于 dW，覆盖余数子断层中心
pymod_m.static_syn(
    output_path="stsyn_ff_zrt.nc",
    src_fault=ff,
    src_fault_size=(1.0, 1.0),
    norths=norths,
    easts=easts,
    calc_upar=True,
)
pymod_m.static_syn(
    output_path="stsyn_ff_zne.nc",
    src_fault=ff,
    src_fault_size=(1.0, 1.0),
    norths=norths,
    easts=easts,
    zne=True,
    calc_upar=True,
)
assert_finite_source_zne("stsyn_ff_zne.nc", "stsyn_ff_zrt.nc")

# The shared Coulomb fixtures cover every supported Kode and the header rake
# marker. Besides checking that the files are accepted, require a finite
# nonzero result so that a silently skipped source cannot pass the test.
coulomb_kodes = Path("cfaults_kodes.inp")
coulomb_rake = Path("cfaults_rake.inr")
coulomb_kodes_suffix = Path("cfaults_kodes_suffix.inr")
coulomb_rake_suffix = Path("cfaults_rake_suffix.inp")
for output_path, src_fault in [
    ("stsyn_ff_kodes.nc", coulomb_kodes),
    ("stsyn_ff_rake.nc", coulomb_rake),
    ("stsyn_ff_kodes_suffix_api.nc", coulomb_kodes_suffix),
    ("stsyn_ff_rake_suffix_api.nc", coulomb_rake_suffix),
]:
    pymod_m.static_syn(
        output_path=output_path,
        src_fault=src_fault,
        src_fault_size=(1.0, 1.0),
        norths=(-4.0, 4.0, 2.0),
        easts=(-4.0, 4.0, 2.0),
    )
    with netcdf_file(output_path, mmap=False) as f:
        z = f.variables["Z"].data
        assert (z == z).all()
        assert abs(z).max() > 0.0

# The header marker, rather than the filename suffix, must select the format.
component, _ = read_static_fields("stsyn_ff_kodes.nc")
component_suffix, _ = read_static_fields("stsyn_ff_kodes_suffix_api.nc")
rake, _ = read_static_fields("stsyn_ff_rake.nc")
rake_suffix, _ = read_static_fields("stsyn_ff_rake_suffix_api.nc")
np.testing.assert_allclose(component_suffix["Z"], component["Z"], rtol=1e-12, atol=1e-12)
np.testing.assert_allclose(rake_suffix["Z"], rake["Z"], rtol=1e-12, atol=1e-12)

# 多台站深度必须给 deprcv
pymod_mr = pygrt.PyModel1D(stgrn="stgrn_mr.nc", modelpath=modname)
pymod_mr.static_greenfn(depsrc=2.0, deprcv=[0.0, 0.5], dists=[0.0, 5.0], calc_upar=True)
pymod_mr.static_syn(scale=1e20, output_path="stsyn_dr.nc", deprcv=0.25, norths=[-2.0, 2.0, 1.0], easts=[-2.0, 2.0, 1.0])

# -------------------- 错误参数 --------------------
try:
    pymod_m.static_syn(scale=1e20, output_path="stsyn_bad.nc", rcv_points=rcv, norths=[-1.0, 1.0, 1.0], easts=[-1.0, 1.0, 1.0])
    raise AssertionError("rcv_points with norths/easts should raise")
except ValueError:
    pass

try:
    pymod_m.static_syn(scale=1e20, output_path="stsyn_bad.nc", rcv_points=rcv, deprcv=0.0)
    raise AssertionError("rcv_points with deprcv should raise")
except ValueError:
    pass

try:
    pymod_m.static_syn(scale=1e20, output_path="stsyn_bad.nc", src_fault=ff)
    raise AssertionError("src_fault with scale should raise")
except ValueError:
    pass

bad_dip = Path("cfaults_bad_dip.inp")
try:
    pymod_m.static_syn(output_path="stsyn_bad.nc", src_fault=bad_dip)
    raise AssertionError("dip=0 should fail in C")
except RuntimeError:
    pass

bad_bot = Path("cfaults_bad_bot.inp")
try:
    pymod_m.static_syn(output_path="stsyn_bad.nc", src_fault=bad_bot)
    raise AssertionError("bot < top should fail in C")
except RuntimeError:
    pass

for name in [
    "stgrn.nc", "stgrn_r.nc", "stgrn_md.nc", "stgrn_mr.nc",
    "stgrn_rf.nc",
    "stgrn_node.nc", "stsyn_node.nc", "stsyn_node_api.nc",
    "stgrn_two_src.nc", "stgrn_two_rcv.nc", "stgrn_two_dist.nc",
    "stsyn_two_src_400.nc", "stsyn_two_src_500.nc", "stsyn_two_rcv_400.nc", "stsyn_two_rcv_500.nc",
    "stsyn_two_dist_400.nc", "stsyn_two_dist_500.nc",
    "stsyn.nc", "stsyn_single_explicit.nc", "stsyn_r.nc",
    "stsyn_md.nc", "stsyn_interp.nc", "stsyn_q.nc", "stsyn_q6.nc", "stsyn_ff.nc",
    "stsyn_rf.nc", "stsyn_rf_default.nc", "stsyn_rq.nc", "stsyn_rf_api.nc",
    "stsyn_ff_zrt.nc", "stsyn_ff_zne.nc", "stsyn_ff_zrt_cli.nc", "stsyn_ff_zne_cli.nc",
    "stsyn_ff_kodes.nc", "stsyn_ff_rake.nc", "stsyn_ff_kodes_suffix.nc", "stsyn_ff_rake_suffix.nc",
    "stsyn_ff_kodes_suffix_api.nc", "stsyn_ff_rake_suffix_api.nc", "stsyn_ff_case.nc",
    "stsyn_dr.nc", "ff_kodes_q.nc",
    "rcv_pts.txt", "rcv_pts_6.txt", "rcv_faults_q.txt", "cfaults_tiny.inp", "cfaults_kodes.inp", "cfaults_rake.inr",
    "cfaults_bad_dip.inp", "cfaults_bad_bot.inp",
]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)

print("test_static_syn.py: all checks passed")
