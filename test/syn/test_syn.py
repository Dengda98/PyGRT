import shutil
from pathlib import Path

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
pymod.compute_grn(
    depsrc=[1.0, depsrc], deprcv=[0.0, deprcv], dists=[5.0, dist],
    nt=nt, dt=dt, calc_upar=True,
)

pymod_root = pymod
pymod = pygrt.PyModel1D(grn="GRN_SINGLE", modelpath=modname)
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, dists=dist, nt=80, dt=dt, calc_upar=True,
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn_single", source="EX",
)
pymod.compute_syn(
    depsrc=depsrc, deprcv=deprcv, dist=dist,
    azimuth=az, scale=1e20, output_path="syn_single_explicit", source="EX",
)

for option, value, message in [
    ("depsrc", depsrc + 0.1, "wrong single source depth should fail"),
    ("deprcv", deprcv + 0.1, "wrong single receiver depth should fail"),
    ("dist", dist + 1.0, "wrong single epicentral distance should fail"),
]:
    try:
        pymod.compute_syn(
            azimuth=az, scale=1e20, output_path="syn_bad", source="EX",
            **{option: value},
        )
        raise AssertionError(message)
    except RuntimeError:
        pass

pymod_root.compute_syn(
    dist=dist, depsrc=depsrc, deprcv=deprcv, azimuth=az, scale=1e20,
    output_path="syn", source="EX",
)
try:
    pygrt.PyModel1D(
        grn="GRN/milrow_2_3_10", modelpath=modname,
    ).compute_syn(
        dist=dist, depsrc=depsrc, deprcv=deprcv, azimuth=az, scale=1e20,
        output_path="syn_subdir_bad", source="EX",
    )
    raise AssertionError("selectors should be rejected for a GF subdirectory")
except RuntimeError:
    pass

# 根目录下存在多个深度和震中距时，-Ds/-Dr/-R 应精确选择目标子目录
for output_path, expected in [
    ("syn_multi_1_0_10", (1.0, 0.0, 10.0)),
    ("syn_multi_2_3_10", (2.0, 3.0, 10.0)),
]:
    sac = read(str(Path(output_path) / "Z.sac"))[0]
    assert abs(sac.stats.sac.evdp - expected[0]) < 1e-5
    assert abs(sac.stats.sac.stel * -1e-3 - expected[1]) < 1e-5
    assert abs(sac.stats.sac.dist - expected[2]) < 1e-5
pymod.compute_syn(
    azimuth=az, scale=1e16, output_path="syn", source="SF",
    force=(-1, 2, -4),
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="DC",
    strike=33, dip=44, rake=55,
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="TS",
    strike=33, dip=44,
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="MT",
    moment_tensor=(1, -2, -5, 0.5, 3, 1.2),
)

# 时间函数
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="EX",
    time_function="p/0.6",
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="EX",
    time_function="t/0.2/0.4/0.7",
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="EX",
    time_function="t/0.4/0.4/0.8",
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="EX",
    time_function="r/1.2",
)

# 积分 / 微分
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="EX",
    integrate_order=1,
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="EX",
    differentiate_order=1,
)

# ZNE / 空间导数
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="EX", zne=True,
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="EX",
    calc_upar=True,
)
pymod.compute_syn(
    azimuth=az, scale=1e20, output_path="syn", source="EX",
    zne=True, calc_upar=True,
)

for name in [
    "GRN", "GRN_SINGLE", "syn", "syn_single", "syn_single_explicit", "syn_bad",
    "syn_subdir_bad",
    "syn_multi_1_0_10", "syn_multi_2_3_10",
]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
