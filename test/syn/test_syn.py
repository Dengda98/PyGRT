import shutil
from pathlib import Path

import numpy as np
from obspy import read

import pygrt

dist = 10.0
depsrc = 2.0
deprcv = 3.0
nt = 600
dt = 0.02
modname = "../milrow"
az = 22.0

pymod = pygrt.PyModel1D(grn="GRN", modelpath=modname)
pymod.greenfn(depsrc=[1.0, depsrc], deprcv=[0.0, deprcv], dists=[5.0, dist], nt=nt, dt=dt, calc_upar=True)

pymod_root = pymod
pymod = pygrt.PyModel1D(grn="GRN_SINGLE", modelpath=modname)
pymod.greenfn(depsrc=depsrc, deprcv=deprcv, dists=dist, nt=80, dt=dt, calc_upar=True)
pymod.syn(azimuth=az, scale=1e20, output_path="syn_single")
pymod.syn(depsrc=depsrc, deprcv=deprcv, dist=dist, azimuth=az, scale=1e20, output_path="syn_single_explicit")

for option, value, message in [
    ("depsrc", depsrc + 0.1, "wrong single source depth should fail"),
    ("deprcv", deprcv + 0.1, "wrong single receiver depth should fail"),
    ("dist", dist + 1.0, "wrong single epicentral distance should fail"),
]:
    try:
        pymod.syn(azimuth=az, scale=1e20, output_path="syn_bad", **{option: value})
        raise AssertionError(message)
    except RuntimeError:
        pass

pymod_root.syn(dist=dist, depsrc=depsrc, deprcv=deprcv, azimuth=az, scale=1e20, output_path="syn")
try:
    pygrt.PyModel1D(grn="GRN/milrow_2_3_10", modelpath=modname).syn(
        dist=dist, depsrc=depsrc, deprcv=deprcv, azimuth=az, scale=1e20,
        output_path="syn_subdir_bad",
    )
    raise AssertionError("selectors should be rejected for a GF subdirectory")
except RuntimeError:
    pass

# 根目录查询位于节点上时，输出几何仍应对应指定位置
for output_path, expected in [
    ("syn_multi_1_0_10", (1.0, 0.0, 10.0)),
    ("syn_multi_2_3_10", (2.0, 3.0, 10.0)),
]:
    sac = read(str(Path(output_path) / "Z.sac"))[0]
    assert abs(sac.stats.sac.evdp - expected[0]) < 1e-5
    assert abs(sac.stats.sac.stel * -1e-3 - expected[1]) < 1e-5
    assert abs(sac.stats.sac.dist - expected[2]) < 1e-5

# 根目录中未落在节点上的单点位置默认插值，也可选择最近邻
for interpolate, output_path in [(True, "syn_root_linear"), (False, "syn_root_nearest")]:
    reference = {path.name: read(str(path))[0].data for path in Path(output_path).glob("*.sac")}
    options = {} if interpolate else dict(interpolate=False)
    pymod_root.syn(depsrc=2, deprcv=3, dist=9, azimuth=az, scale=1e20, output_path=output_path, **options)
    for name, data in reference.items():
        np.testing.assert_array_equal(read(str(Path(output_path)/name))[0].data, data)
assert not np.array_equal(read("syn_root_linear/Z.sac")[0].data, read("syn_root_nearest/Z.sac")[0].data)

# 节点子目录的插值选项不起作用，与同位置的根目录合成一致
node_model = pygrt.PyModel1D(grn="GRN/milrow_2_3_10/.", modelpath=modname)
reference = read("syn_multi_2_3_10/Z.sac")[0].data
for interpolate in [False, True]:
    output_path = f"syn_node_i{int(interpolate)}"
    np.testing.assert_array_equal(read(f"{output_path}/Z.sac")[0].data, reference)
    node_model.syn(azimuth=az, scale=1e20, interpolate=interpolate, output_path=output_path)
    np.testing.assert_array_equal(read(f"{output_path}/Z.sac")[0].data, reference)
for option, value in [("depsrc", 2), ("deprcv", 3), ("dist", 10)]:
    try:
        node_model.syn(azimuth=az, scale=1e20, output_path="syn_subdir_bad", **{option: value})
    except RuntimeError:
        pass
    else:
        raise AssertionError(f"A node directory must reject {option}")
pymod.syn(azimuth=az, scale=1e16, output_path="syn", force=(-1, 2, -4))
pymod.syn(azimuth=az, scale=1e20, output_path="syn", strike=33, dip=44, rake=55)
pymod.syn(azimuth=az, scale=1e20, output_path="syn", strike=33, dip=44)
pymod.syn(azimuth=az, scale=1e20, output_path="syn", moment_tensor=(1, -2, -5, 0.5, 3, 1.2))

# 震源时间函数按矩形法进行面积归一化
# 自定义时间函数检查 dt 乘样本和，非单位面积时警告并自动归一化


custom_signal = read("syn_custom/sig.sac")[0].data
assert np.isclose(np.sum(custom_signal), 1.0 / dt, rtol=1e-5, atol=1e-5)
custom_warning_signal = read("syn_custom_warning/sig.sac")[0].data
assert np.min(custom_warning_signal) < 0.0
assert np.isclose(np.sum(custom_warning_signal), 1.0 / dt, rtol=1e-5, atol=1e-5)

for time_function in ["p/0.6", "t/0.2/0.2/0.3", "t/0.4/0/0.4", "c/0.2/0.4"]:
    pymod.syn(azimuth=az, scale=1e20, output_path="syn", time_function=time_function)
    trace = read("syn/sig.sac")[0]
    assert np.isclose(np.sum(trace.data)*dt, 1.0, rtol=1e-5, atol=1e-5)

pymod.syn(azimuth=az, scale=1e20, output_path="syn", time_function="r/1.2")
signal = read("syn/sig.sac")[0].data
assert np.max(signal) <= 1.0 + 1e-5

# 积分 / 微分
pymod.syn(azimuth=az, scale=1e20, output_path="syn", integrate_order=1)
pymod.syn(azimuth=az, scale=1e20, output_path="syn", differentiate_order=1)

# ZNE / 空间导数
pymod.syn(azimuth=az, scale=1e20, output_path="syn", zne=True)
pymod.syn(azimuth=az, scale=1e20, output_path="syn", calc_upar=True)
pymod.syn(azimuth=az, scale=1e20, output_path="syn", zne=True, calc_upar=True)

# 单台动态应力投影在水平面上应退化为 ZZ 法向应力和 ZN 剪应力
pygrt.utils.stress("syn")
pygrt.utils.sproj("syn", strike=0, dip=0, rake=0)
normal = read("syn/sigma_n.sac")[0].data
shear = read("syn/tau_s.sac")[0].data
np.testing.assert_array_equal(normal, read("syn/stress_ZZ.sac")[0].data)
np.testing.assert_array_equal(shear, read("syn/stress_ZN.sac")[0].data)
pygrt.utils.coulomb("syn", 0.4)
np.testing.assert_allclose(read("syn/coulomb.sac")[0].data, shear + 0.4*normal, rtol=2e-6, atol=np.max(np.abs(normal))*1e-7)

# 同一接收面在 ZRT 和 ZNE 输入下应得到一致的投影
pygrt.utils.sproj("syn", strike=33, dip=44, rake=55)
reference = {name: read(f"syn/{name}.sac")[0].data for name in ["sigma_n", "tau_s"]}
pymod.syn(azimuth=az, scale=1e20, output_path="syn_projection_zrt", calc_upar=True)
pygrt.utils.stress("syn_projection_zrt")
pygrt.utils.sproj("syn_projection_zrt", strike=33, dip=44, rake=55)
for name, expected in reference.items():
    np.testing.assert_allclose(read(f"syn_projection_zrt/{name}.sac")[0].data, expected,
                               rtol=2e-5, atol=np.max(np.abs(expected))*2e-7)
shutil.rmtree("syn_projection_zrt")

for name in [
    "GRN", "GRN_SINGLE", "syn", "syn_single", "syn_single_explicit", "syn_bad",
    "syn_custom", "syn_custom_warning",
    "syn_subdir_bad",
    "syn_multi_1_0_10", "syn_multi_2_3_10",
    "syn_root_linear", "syn_root_nearest", "syn_node_i0", "syn_node_i1",
]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
