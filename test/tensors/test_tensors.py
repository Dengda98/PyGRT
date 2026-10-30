import shutil
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

model = pygrt.PyModel1D(grn="GRN", modelpath="../milrow")
model.greenfn(depsrc=2, deprcv=[0, 1], dists=[1, 6, 12], nt=16, dt=0.1, calc_upar=True, print_log=False)
model.syn(depsrc=2, deprcv=0, dist=6, azimuth=22, scale=1e20, output_path="syn_plain")
# 未计算位移导数，不可计算应变
with raises(RuntimeError):
    pygrt.tensor.strain("syn_plain")

model.syn(depsrc=2, deprcv=0, dist=6, azimuth=22, scale=1e20, calc_upar=True, output_path="syn")
pygrt.tensor.strain("syn")
pygrt.tensor.stress("syn")
pygrt.tensor.rotation("syn")
# 未指定接收断层几何
with raises(RuntimeError):
    pygrt.tensor.sproj("syn")
# 未计算法向和剪切应力投影，不可计算 Coulomb 应力
with raises(RuntimeError):
    pygrt.tensor.coulomb("syn", 0.4)
pygrt.tensor.sproj("syn", strike=33, dip=44, rake=55)
pygrt.tensor.coulomb("syn", 0.4)

model.syn(depsrc=2, deprcv=0, dist=6, azimuth=22, scale=1e20, zne=True, calc_upar=True, output_path="syn")
pygrt.tensor.strain("syn")
pygrt.tensor.stress("syn")
pygrt.tensor.rotation("syn")
pygrt.tensor.sproj("syn", strike=33, dip=44, rake=55)
pygrt.tensor.coulomb("syn", 0)

model.syn(depsrc=2, scale=1e20, rcv_points="rcv_geometry.txt", calc_upar=True, output_path="syn_points")
pygrt.tensor.strain("syn_points")
pygrt.tensor.stress("syn_points")
pygrt.tensor.rotation("syn_points")
pygrt.tensor.sproj("syn_points")
pygrt.tensor.sproj("syn_points", rcv_points="rcv_geometry.txt")
pygrt.tensor.coulomb("syn_points", 0.4)

model.syn(depsrc=2, scale=1e20, rcv_fault="rcv_faults.inr", calc_upar=True, output_path="syn_faults")
pygrt.tensor.stress("syn_faults")
pygrt.tensor.sproj("syn_faults")
pygrt.tensor.sproj("syn_faults", rake=55, force_rake=True)
pygrt.tensor.coulomb("syn_faults", 0.4)

# 未指定接收断层的滑动角
with raises(ValueError):
    pygrt.tensor.sproj("syn", strike=33, dip=44)
# 摩擦系数不可为负数
with raises(RuntimeError):
    pygrt.tensor.coulomb("syn", -0.1)
# 输入的动态合成目录不存在
with raises(FileNotFoundError):
    pygrt.tensor.stress("missing")

for name in ["GRN", "syn", "syn_points", "syn_faults", "syn_plain"]:
    shutil.rmtree(name)
