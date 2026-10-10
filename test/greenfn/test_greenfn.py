import shutil
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

model = pygrt.PyModel1D(grn="GRN", modelpath="../milrow")
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1)
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1, calc_upar=True, upsampling_n=2, zeta=0.6, nthreads=2)
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1, freqband=[0.5, 2], Length=10, keepAllFreq=True, skipImagComps=True)
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1, k0=4, use_kmax_ref=True, ampk=1.2, keps=1e-3, vmin_ref=1.5)
model.greenfn(depsrc=2, deprcv=0, dists=100, nt=32, dt=1, filonLength=10, filonCut=2, delayT0=-2, delayV0=9)
model.greenfn(depsrc=2, deprcv=0, dists=100, nt=32, dt=1, safilonTol=1e-3, ref_first_p=True, delayT0=-1)
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1, converg_method="DCM")
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1, converg_method="PTAM")
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1, converg_method="none", gf_source=["EX", "VF"], statsidxs=[1, 2])
model.greenfn(depsrc=[1, 2], deprcv=[0, 1], dists=[2, 5, 10], nt=16, dt=0.1, print_log=False)

model = pygrt.PyModel1D(grn="GRN", modelpath="../milrow", topbound="rigid", botbound="free")
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1)
model = pygrt.PyModel1D(grn="GRN", modelpath="../milrow", topbound="halfspace", botbound="rigid")
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1)

# 频带下限不可超过上限
with raises(RuntimeError):
    model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1, freqband=[2, 0.5], print_log=False)
# 不可同时设置 FIM 和 SAFIM
with raises(RuntimeError):
    model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1, filonLength=10, safilonTol=1e-3)

# 震中距列表不可包含逆序或重复值
with raises(ValueError):
    model.greenfn(depsrc=2, deprcv=0, dists=[5, 2], nt=32, dt=0.1)
# 时间采样间隔不可为零或负数
with raises(ValueError):
    model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0)
# 不可设置不支持的收敛方法
with raises(ValueError):
    model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=0.1, converg_method="bad")

shutil.rmtree("GRN")
shutil.rmtree("GRN_grtstats")
