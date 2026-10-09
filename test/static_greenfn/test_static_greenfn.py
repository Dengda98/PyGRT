import shutil
from pathlib import Path
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

model = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="../milrow")
model.static_greenfn(depsrc=2, deprcv=0, dists=[0, 2, 4])
model.static_greenfn(depsrc=2, deprcv=0, dists=[0, 2, 4], calc_upar=True, Length=10)
model.static_greenfn(depsrc=2, deprcv=0, dists=[0, 2, 4], k0=4, use_kmax_ref=True, keps=1e-3, stats=True)
model.static_greenfn(depsrc=2, deprcv=0, dists=[0, 2, 4], converg_method="DCM")
model.static_greenfn(depsrc=2, deprcv=0, dists=[0, 2, 4], converg_method="PTAM")
model.static_greenfn(depsrc=2, deprcv=0, dists=[0, 2, 4], converg_method="none")
model.static_greenfn(depsrc=2, deprcv=0, dists=[100, 120], filonLength=10, filonCut=2)
model.static_greenfn(depsrc=2, deprcv=0, dists=[100, 120], safilonTol=1e-3)
model.static_greenfn(depsrc=2, deprcv=0, norths=[-2, 2, 2], easts=[-2, 2, 2])
model.static_greenfn(depsrc=[1, 2], deprcv=[0, 1], dists=[2, 5, 10], calc_upar=True)

model = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="../milrow", topbound="rigid", botbound="free")
model.static_greenfn(depsrc=2, deprcv=0, dists=[0, 2, 4])
model = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="../milrow", topbound="halfspace", botbound="rigid")
model.static_greenfn(depsrc=2, deprcv=0, dists=[0, 2, 4])

# 接收网格不可只设置北向坐标
with raises(ValueError):
    model.static_greenfn(depsrc=2, deprcv=0, norths=[-2, 2, 2])

# 震中距列表不可包含逆序或重复值
with raises(ValueError):
    model.static_greenfn(depsrc=2, deprcv=0, dists=[5, 2])
# 不可同时设置震中距和接收网格
with raises(ValueError):
    model.static_greenfn(depsrc=2, deprcv=0, dists=5, norths=[0, 2, 1], easts=[0, 2, 1])
# 震源深度不可为负数
with raises(ValueError):
    model.static_greenfn(depsrc=-1, deprcv=0, dists=5)

Path("stgrn.nc").unlink()
shutil.rmtree("stgrtstats")
