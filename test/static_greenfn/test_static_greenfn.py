import shutil
import warnings
from pathlib import Path

import pygrt
from scipy.io import netcdf_file

depsrc = 2.0
deprcv = 0.0
# norths/easts 各为 start/stop/step (km)
norths = [-3.0, 3.0, 0.2]
easts = [-2.0, 2.0, 0.2]
modname = "../milrow"

pymod = pygrt.PyModel1D(modname)
pymod.set_static_grn_path("stgrn.nc")

pymod.compute_static_grn(depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, calc_upar=True,
)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, Length=20,
)

pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts,
    Length=20, converg_method="DCM",
)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts,
    Length=20, converg_method="PTAM",
)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts,
    Length=20, converg_method="none",
)

pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, k0=4, keps=1e-3,
)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, k0=4, stats=True,
)

# 单深度输出应为 4D STGRNLIB（各深度维长度为 1）
with netcdf_file("stgrn.nc", mmap=False) as f:
    assert "depsrc" in f.dimensions and "deprcv" in f.dimensions
    assert f.dimensions["depsrc"] == 1
    assert f.dimensions["deprcv"] == 1

# boundary condition
pymod = pygrt.PyModel1D(modname, topbound="free", botbound="free")
pymod.set_static_grn_path("stgrn.nc")
pymod.compute_static_grn(depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)

pymod = pygrt.PyModel1D(modname, topbound="halfspace", botbound="free")
pymod.set_static_grn_path("stgrn.nc")
pymod.compute_static_grn(depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)

pymod = pygrt.PyModel1D(modname, topbound="rigid", botbound="rigid")
pymod.set_static_grn_path("stgrn.nc")
pymod.compute_static_grn(depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)

# -------------------- 多深度功能 --------------------
norths_c = [-2.0, 2.0, 1.0]
easts_c = [-2.0, 2.0, 1.0]
depsrcs = [1.0, 2.0, 3.0]
deprcvs = [0.0, 0.5]

pymod_m = pygrt.PyModel1D(modname)
pymod_m.set_static_grn_path("stgrn_py_multi.nc")
pymod_m.compute_static_grn(depsrc=depsrcs, deprcv=deprcvs, norths=norths_c, easts=easts_c)
assert Path("stgrn_py_multi.nc").is_file()
with netcdf_file("stgrn_py_multi.nc", mmap=False) as f:
    assert f.dimensions["depsrc"] == 3
    assert f.dimensions["deprcv"] == 2

# 仅多震源深度
pymod_ms = pygrt.PyModel1D(modname)
pymod_ms.set_static_grn_path("stgrn_py_ms.nc")
pymod_ms.compute_static_grn(
    depsrc=[1.0, 2.0], deprcv=0.0, norths=norths_c, easts=easts_c, calc_upar=True,
)

# 仅多台站深度
pymod_mr = pygrt.PyModel1D(modname)
pymod_mr.set_static_grn_path("stgrn_py_mr.nc")
pymod_mr.compute_static_grn(depsrc=2.0, deprcv=[0.0, 0.5], norths=norths_c, easts=easts_c)

# -R / distarr 建库
pymod_r = pygrt.PyModel1D(modname)
pymod_r.set_static_grn_path("stgrn_py_r.nc")
pymod_r.compute_static_grn(
    depsrc=2.0, deprcv=0.0, distarr=[0.0, 1.0, 2.0, 4.0],
)
with netcdf_file("stgrn_py_r.nc", mmap=False) as f:
    assert f.dimensions["north"] == 1
    assert f.dimensions["east"] == 4

# -------------------- 错误 / 警告（Python）--------------------
try:
    pymod_m.compute_static_grn(depsrc=-1.0, deprcv=0.0, norths=norths_c, easts=easts_c)
    raise AssertionError("negative depsrc should raise")
except ValueError:
    pass

try:
    pymod_m.compute_static_grn(depsrc=1.0, deprcv=-0.5, norths=norths_c, easts=easts_c)
    raise AssertionError("negative deprcv should raise")
except ValueError:
    pass

try:
    pymod_m.compute_static_grn(depsrc=[], deprcv=0.0, norths=norths_c, easts=easts_c)
    raise AssertionError("empty depsrc should raise")
except ValueError:
    pass

try:
    pymod_r.compute_static_grn(depsrc=2.0, deprcv=0.0, distarr=[0.0, 2.0, 1.0])
    raise AssertionError("non-ascending distarr should raise")
except ValueError:
    pass

# 多深度 stats 应警告并忽略
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    pymod_m.set_static_grn_path("stgrn_py_multi2.nc")
    pymod_m.compute_static_grn(
        depsrc=depsrcs, deprcv=deprcvs, norths=norths_c, easts=easts_c, stats=True,
    )
    assert any("stats" in str(x.message) for x in w)

for name in [
    "stgrn.nc", "stgrn_py_multi.nc", "stgrn_py_ms.nc", "stgrn_py_mr.nc",
    "stgrn_py_multi2.nc", "stgrn_py_r.nc", "stgrtstats",
]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)

print("test_static_greenfn.py: all checks passed")
