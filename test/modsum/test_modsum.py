import shutil
from pathlib import Path
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

model = pygrt.PyModel1D(grn="GRN", modelpath="../milrow")
model.eigenv(wtype="R", freqs=(0, 0.5, 0.1), phase_path="phase_R.nc", all_modes=True)
model.eigenv(wtype="L", freqs=(0, 0.5, 0.1), phase_path="phase_L.nc", all_modes=True)
model.modsum(phase_path="phase_R.nc", depsrc=2, deprcv=1, dists=100, modes=0, upsampling_n=2, calc_upar=True)
model.modsum(phase_path="phase_L.nc", depsrc=2, deprcv=1, dists=100, all_modes=True, output_path="GRN", upsampling_n=2)
model.modsum(phase_path="phase_R.nc", depsrc=[1, 2], deprcv=[0, 1], dists=[100, 120], modes=0,
             freqband=(0.1, 0.4), delayT0=-2, delayV0=9, nthreads=2)
model.modsum(phase_path="phase_R.nc", depsrc=2, deprcv=1, dists=100, modes=0, ref_first_p=True, delayT0=-1, gf_source=["EX", "VF"])

# 升采样倍数不可为零或负数
with raises(ValueError):
    model.modsum(phase_path="phase_R.nc", depsrc=2, deprcv=1, dists=100, modes=0, upsampling_n=0)

# 不可同时设置特定阶数和全部阶数
with raises(ValueError):
    model.modsum(phase_path="phase_R.nc", depsrc=2, deprcv=1, dists=100, modes=0, all_modes=True)
# 震中距列表不可包含逆序或重复值
with raises(ValueError):
    model.modsum(phase_path="phase_R.nc", depsrc=2, deprcv=1, dists=[100, 50], modes=0)

shutil.rmtree("GRN")
Path("phase_R.nc").unlink()
Path("phase_L.nc").unlink()
