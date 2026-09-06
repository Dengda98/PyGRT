from contextlib import redirect_stdout
from pathlib import Path

import pygrt

modname = "../milrow"
pymod = pygrt.PyModel1D(modelpath=modname)

# Rayleigh / Love 全阶频散
pymod.eigenv(wtype="R", freqs=(0.0, 0.5, 0.05), phase_path="phase_R.nc", all_modes=True)
pymod.eigenv(wtype="L", freqs=(0.0, 0.5, 0.05), phase_path="phase_L.nc", all_modes=True)
assert Path("phase_R.nc").is_file()
assert Path("phase_L.nc").is_file()

# 单频、周期形式、最大阶数、搜根控制参数
pymod.eigenv(wtype="R", freqs=1.0, phase_path="phase_R1.nc")
pymod.eigenv(wtype="R", periods=(1.0, 5.0, 1.0), phase_path="phase_Rp.nc", max_order=2)
pymod.eigenv(wtype="L", freqs=(0.0, 0.5, 0.05), phase_path="phase_Lt.nc", all_modes=True, ctrl_kw={"tol": 3.0})

# 久期函数 debug，输出重定向到文件
with open("secfunc_R.txt", "w") as f, redirect_stdout(f):
    pymod.eigenv(wtype="R", secular_freq=1.0, cmin=3.0, cmax=5.5, iref=0)
assert Path("secfunc_R.txt").stat().st_size > 0

try:
    pymod.eigenv(wtype="X", freqs=1.0, phase_path="bad.nc")
except ValueError:
    pass
else:
    raise AssertionError("invalid wtype should raise")

try:
    pymod.eigenv(wtype="R", freqs=1.0, periods=2.0, phase_path="bad.nc")
except ValueError:
    pass
else:
    raise AssertionError("freqs and periods should be mutually exclusive")

try:
    pymod.eigenv(wtype="R", freqs=1.0, phase_path="bad.nc", secular_freq=1.0)
except ValueError:
    pass
else:
    raise AssertionError("secular_freq mixed with freqs should raise")

try:
    pymod.eigenv(wtype="R", freqs=1.0, phase_path="bad.nc", ctrl_kw={"foo": 1.0})
except ValueError:
    pass
else:
    raise AssertionError("unknown ctrl_kw key should raise")

for name in ["phase_R.nc", "phase_L.nc", "phase_R1.nc", "phase_Rp.nc", "phase_Lt.nc", "secfunc_R.txt"]:
    Path(name).unlink(missing_ok=True)
