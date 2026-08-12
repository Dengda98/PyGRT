import shutil
from pathlib import Path

import pygrt

dist = 10.0
depsrc = 2.0
deprcv = 3.0
nt = 600
dt = 0.02
modname = "../milrow"
az = 22.0

pymod = pygrt.PyModel1D(modname)
pymod.set_dynamic_grn_path("GRN")
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt, calc_upar=True,
)

pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e16, output_path="syn", source="SF",
    force=(-1, 2, -4),
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="DC",
    strike=33, dip=44, rake=55,
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="TS",
    strike=33, dip=44,
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="MT",
    moment_tensor=(1, -2, -5, 0.5, 3, 1.2),
)

# 时间函数
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
    time_function="p/0.6",
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
    time_function="t/0.2/0.4/0.7",
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
    time_function="t/0.4/0.4/0.8",
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
    time_function="r/1.2",
)

# 积分 / 微分
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
    integrate_order=1,
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
    differentiate_order=1,
)

# ZNE / 空间导数
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX", zne=True,
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
    calc_upar=True,
)
pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
    zne=True, calc_upar=True,
)

for name in ["GRN", "syn"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
