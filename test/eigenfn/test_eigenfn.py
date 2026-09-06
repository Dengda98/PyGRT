from pathlib import Path

import pygrt

pymod = pygrt.PyModel1D(modelpath="../milrow")
pymod.eigenv(wtype="R", freqs=(0.0, 0.5, 0.05), phase_path="phase_R.nc", all_modes=True)

# 本征函数、群速度、敏感核、能量积分
pymod.eigenfn(
    phase_path="phase_R.nc",
    freqs=0.5,
    modes=(0, 10),
    eigenfn_path="egn_R.nc",
    depths=(0.0, 20.0, 0.5),
    group_path="group_R.nc",
    csens_path="csens.nc",
    usens_path="usens.nc",
    egy_path="egyint.nc",
    sensitivity_dz=0.5,
)
assert Path("egn_R.nc").is_file()
assert Path("group_R.nc").is_file()
assert Path("csens.nc").is_file()
assert Path("usens.nc").is_file()
assert Path("egyint.nc").is_file()

# 仅输出群速度
pymod.eigenfn(phase_path="phase_R.nc", all_modes=True, group_path="group_all.nc")
assert Path("group_all.nc").is_file()

# 周期选择
pymod.eigenfn(phase_path="phase_R.nc", periods=2.0, modes=0, eigenfn_path="egn_p.nc", depths=1.0)
assert Path("egn_p.nc").is_file()

# 频率范围选择
pymod.eigenfn(phase_path="phase_R.nc", freqs=(0.1, 0.4), modes=0, group_path="group_band.nc")
assert Path("group_band.nc").is_file()

try:
    pymod.eigenfn(phase_path="phase_R.nc", modes=0, all_modes=True)
except ValueError:
    pass
else:
    raise AssertionError("modes and all_modes should be mutually exclusive")

try:
    pymod.eigenfn(phase_path="phase_R.nc", freqs=1.0, periods=2.0)
except ValueError:
    pass
else:
    raise AssertionError("freqs and periods should be mutually exclusive")

try:
    pymod.eigenfn(phase_path="phase_R.nc", depths=0.0)
except ValueError:
    pass
else:
    raise AssertionError("depths without eigenfn_path should raise")

for name in ["phase_R.nc", "egn_R.nc", "group_R.nc", "csens.nc", "usens.nc", "egyint.nc", "group_all.nc", "egn_p.nc", "group_band.nc"]:
    Path(name).unlink(missing_ok=True)
