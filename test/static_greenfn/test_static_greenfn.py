import shutil
import warnings
from pathlib import Path

import pygrt
from scipy.io import netcdf_file

depsrc = 2.0
deprcv = 0.0
dists = [0.0, 1.0, 2.0, 3.0, 4.0]
modname = "../milrow"

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=modname)

pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists)
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists, calc_upar=True)
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists, Length=20)

pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists, Length=20, converg_method="DCM")
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists, Length=20, converg_method="PTAM")
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists, Length=20, converg_method="none")

pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists, k0=4, keps=1e-3)
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists, k0=4, stats=True)

# -X/-Y 仅测试二维网格输入
pymod.static_greenfn(
    depsrc=depsrc,
    deprcv=deprcv,
    norths=[-2.0, 2.0, 1.0],
    easts=[-2.0, 2.0, 1.0],
)

# 单深度输出应为 4D STGRNLIB（各深度维长度为 1）
with netcdf_file("stgrn.nc", mmap=False) as f:
    assert "depsrc" in f.dimensions and "deprcv" in f.dimensions
    assert f.dimensions["depsrc"] == 1
    assert f.dimensions["deprcv"] == 1

# boundary condition
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=modname, topbound="free", botbound="free")
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists)

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=modname, topbound="halfspace", botbound="free")
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists)

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=modname, topbound="rigid", botbound="rigid")
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, dists=dists)

# -------------------- 多深度功能 --------------------
depsrcs = [1.0, 2.0, 3.0]
deprcvs = [0.0, 0.5]

pymod_m = pygrt.PyModel1D(stgrn="stgrn_py_multi.nc", modelpath=modname)
pymod_m.static_greenfn(depsrc=depsrcs, deprcv=deprcvs, dists=dists)
assert Path("stgrn_py_multi.nc").is_file()
with netcdf_file("stgrn_py_multi.nc", mmap=False) as f:
    assert f.dimensions["depsrc"] == 3
    assert f.dimensions["deprcv"] == 2

# 仅多震源深度
pymod_ms = pygrt.PyModel1D(stgrn="stgrn_py_ms.nc", modelpath=modname)
pymod_ms.static_greenfn(depsrc=[1.0, 2.0], deprcv=0.0, dists=dists, calc_upar=True)

# 仅多台站深度
pymod_mr = pygrt.PyModel1D(stgrn="stgrn_py_mr.nc", modelpath=modname)
pymod_mr.static_greenfn(depsrc=2.0, deprcv=[0.0, 0.5], dists=dists)

# -R / dists 建库
pymod_r = pygrt.PyModel1D(stgrn="stgrn_py_r.nc", modelpath=modname)
pymod_r.static_greenfn(depsrc=2.0, deprcv=0.0, dists=[0.0, 1.0, 2.0, 4.0])
with netcdf_file("stgrn_py_r.nc", mmap=False) as f:
    assert f.dimensions["north"] == 1
    assert f.dimensions["east"] == 4

# -------------------- 错误 / 警告（Python）--------------------
try:
    pymod_m.static_greenfn(depsrc=-1.0, deprcv=0.0, dists=dists)
    raise AssertionError("negative depsrc should raise")
except ValueError:
    pass

try:
    pymod_m.static_greenfn(depsrc=1.0, deprcv=-0.5, dists=dists)
    raise AssertionError("negative deprcv should raise")
except ValueError:
    pass

try:
    pymod_m.static_greenfn(depsrc=[], deprcv=0.0, dists=dists)
    raise AssertionError("empty depsrc should raise")
except ValueError:
    pass

try:
    pymod_r.static_greenfn(depsrc=2.0, deprcv=0.0, dists=[0.0, 2.0, 1.0])
    raise AssertionError("non-ascending dists should raise")
except ValueError:
    pass

# 多深度 stats 应警告并忽略
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    pymod_m = pygrt.PyModel1D(stgrn="stgrn_py_multi2.nc", modelpath=modname)
    pymod_m.static_greenfn(depsrc=depsrcs, deprcv=deprcvs, dists=dists, stats=True)
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
