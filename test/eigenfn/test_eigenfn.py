from pathlib import Path
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

model = pygrt.PyModel1D(modelpath="../milrow")
model.eigenv(wtype="R", freqs=(0, 0.5, 0.1), phase_path="phase.nc", all_modes=True)
model.eigenfn(phase_path="phase.nc", freqs=0.5, modes=(0, 2), eigenfn_path="egn.nc", depths=(0, 10, 1),
              group_path="group.nc", csens_path="csens.nc", usens_path="usens.nc", egy_path="energy.nc", sensitivity_dz=1)
model.eigenfn(phase_path="phase.nc", periods=2, modes=0, eigenfn_path="egn.nc", depths=1)
model.eigenfn(phase_path="phase.nc", freqs=(0.2, 0.5), all_modes=True, group_path="group.nc")

# 模态阶数不可为负数
with raises(ValueError):
    model.eigenfn(phase_path="phase.nc", modes=-1, group_path="group.nc")

# 不可同时设置特定阶数和全部阶数
with raises(ValueError):
    model.eigenfn(phase_path="phase.nc", modes=0, all_modes=True)
# 不可同时设置频率和周期
with raises(ValueError):
    model.eigenfn(phase_path="phase.nc", freqs=0.5, periods=2)
# 未指定本征函数输出文件，不可设置深度
with raises(ValueError):
    model.eigenfn(phase_path="phase.nc", depths=1)

for name in ["phase.nc", "egn.nc", "group.nc", "csens.nc", "usens.nc", "energy.nc"]:
    Path(name).unlink()
