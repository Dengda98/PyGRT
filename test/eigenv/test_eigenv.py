from pathlib import Path
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

model = pygrt.PyModel1D(modelpath="../milrow")
model.eigenv(wtype="R", freqs=(0, 0.5, 0.1), phase_path="phase_R.nc", all_modes=True)
model.eigenv(wtype="L", freqs=(0, 0.5, 0.1), phase_path="phase_L.nc", max_order=2, ctrl_kw={"tol": 3}, nthreads=2)
model.eigenv(wtype="R", periods=(2, 4, 1), phase_path="phase_p.nc")
model.eigenv(wtype="R", secular_freq=1, cmin=3, cmax=5.5, iref=0)

# 周期不可为零或负数
with raises(ValueError):
    model.eigenv(wtype="R", periods=0, phase_path="phase_bad.nc")

# 面波类型不可设置为 R、L 之外的值
with raises(ValueError):
    model.eigenv(wtype="X", freqs=0.5, phase_path="phase_bad.nc")
# 不可同时设置频率和周期
with raises(ValueError):
    model.eigenv(wtype="R", freqs=0.5, periods=2, phase_path="phase_bad.nc")
# 不可设置不支持的搜根控制参数
with raises(ValueError):
    model.eigenv(wtype="R", freqs=0.5, phase_path="phase_bad.nc", ctrl_kw={"bad": 1})

for path in Path(".").glob("phase_*.nc"):
    path.unlink()
