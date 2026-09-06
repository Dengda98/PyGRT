import shutil
from pathlib import Path

import pygrt

pymod = pygrt.PyModel1D(grn="GRN_NM_0", modelpath="../milrow")
pymod.eigenv(wtype="R", freqs=(0.0, 0.5, 0.05), phase_path="phase_R.nc", all_modes=True)
pymod.eigenv(wtype="L", freqs=(0.0, 0.5, 0.05), phase_path="phase_L.nc", all_modes=True)

# 未指定 output_path 时写入构造时的 grn 根目录
pymod.modsum(phase_path="phase_R.nc", depsrc=2.0, deprcv=1.0, dists=100.0, modes=0, upsampling_n=5, calc_upar=True)
assert Path("GRN_NM_0/milrow").is_file()
pymod.modsum(phase_path="phase_L.nc", depsrc=2.0, deprcv=1.0, dists=100.0, modes=0, upsampling_n=5, calc_upar=True)

# 指定输出目录、其它阶数
pymod.modsum(phase_path="phase_R.nc", depsrc=2.0, deprcv=1.0, dists=100.0, modes=1, output_path="GRN_NM_1", upsampling_n=5)
pymod.modsum(phase_path="phase_R.nc", depsrc=2.0, deprcv=1.0, dists=100.0, all_modes=True, output_path="GRN_NM_all", upsampling_n=5)

# 多震源/台站深度
pymod.modsum(
    phase_path="phase_R.nc",
    depsrc=[1.0, 2.0],
    deprcv=[0.0, 1.0],
    dists=100.0,
    modes=0,
    output_path="GRN_NM_MULTI",
    upsampling_n=2,
    calc_upar=True,
)
assert Path("GRN_NM_MULTI/milrow_1_0_100/EXZ.sac").is_file()
assert Path("GRN_NM_MULTI/milrow_2_1_100/EXZ.sac").is_file()

# 频带截取
pymod.modsum(phase_path="phase_R.nc", depsrc=2.0, deprcv=1.0, dists=100.0, modes=0, output_path="GRN_NM_F", freqband=(0.1, 0.4))

try:
    pymod.modsum(phase_path="phase_R.nc", depsrc=2.0, deprcv=1.0, dists=100.0, modes=0, all_modes=True)
except ValueError:
    pass
else:
    raise AssertionError("modes and all_modes should be mutually exclusive")

for name in ["GRN_NM_0", "GRN_NM_1", "GRN_NM_all", "GRN_NM_MULTI", "GRN_NM_F", "phase_R.nc", "phase_L.nc"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
