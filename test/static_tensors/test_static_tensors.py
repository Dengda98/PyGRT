import subprocess
from pathlib import Path
from unittest import TestCase

import pygrt

subprocess.run(["bash", "-ec", "source ../common.sh; create_test_files"], check=True)

raises = TestCase().assertRaises

model = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="../milrow")
model.static_greenfn(depsrc=2, deprcv=[0, 1], dists=[0, 3, 6, 9, 12], calc_upar=True)
model.static_syn(scale=1e20, deprcv=0, output_path="stsyn_plain.nc")
# 未计算位移导数，不可计算应变
with raises(RuntimeError):
    pygrt.utils.static_strain("stsyn_plain.nc")

model.static_syn(scale=1e20, deprcv=0, norths=[-2, 2, 2], easts=[-2, 2, 2], calc_upar=True, output_path="stsyn.nc")
pygrt.utils.static_strain("stsyn.nc")
pygrt.utils.static_stress("stsyn.nc")
pygrt.utils.static_rotation("stsyn.nc")
# 未指定接收断层几何
with raises(RuntimeError):
    pygrt.utils.static_sproj("stsyn.nc")
# 未计算法向和剪切应力投影，不可计算 Coulomb 应力
with raises(RuntimeError):
    pygrt.utils.static_coulomb("stsyn.nc", 0.4)
pygrt.utils.static_sproj("stsyn.nc", strike=33, dip=44, rake=55)
pygrt.utils.static_coulomb("stsyn.nc", 0.4)

model.static_syn(scale=1e20, deprcv=0, zne=True, calc_upar=True, output_path="stsyn.nc")
pygrt.utils.static_strain("stsyn.nc")
pygrt.utils.static_stress("stsyn.nc")
pygrt.utils.static_rotation("stsyn.nc")
pygrt.utils.static_sproj("stsyn.nc", strike=33, dip=44, rake=55)
pygrt.utils.static_coulomb("stsyn.nc", 0)

model.static_syn(scale=1e20, rcv_points="rcv_geometry.txt", calc_upar=True, output_path="stsyn_points.nc")
pygrt.utils.static_strain("stsyn_points.nc")
pygrt.utils.static_stress("stsyn_points.nc")
pygrt.utils.static_rotation("stsyn_points.nc")
pygrt.utils.static_sproj("stsyn_points.nc")
pygrt.utils.static_sproj("stsyn_points.nc", rcv_points="rcv_geometry.txt")
pygrt.utils.static_coulomb("stsyn_points.nc", 0.4)

model.static_syn(scale=1e20, rcv_fault="rcv_faults.inr", calc_upar=True, output_path="stsyn_faults.nc")
pygrt.utils.static_stress("stsyn_faults.nc")
pygrt.utils.static_sproj("stsyn_faults.nc")
pygrt.utils.static_sproj("stsyn_faults.nc", rake=55, force_rake=True)
pygrt.utils.static_coulomb("stsyn_faults.nc", 0.4)

# 未指定接收断层的滑动角
with raises(ValueError):
    pygrt.utils.static_sproj("stsyn.nc", strike=33, dip=44)
# 摩擦系数不可为负数
with raises(RuntimeError):
    pygrt.utils.static_coulomb("stsyn.nc", -0.1)
# 输入的静态合成文件不存在
with raises(FileNotFoundError):
    pygrt.utils.static_strain("missing.nc")

for name in ["stgrn.nc", "stsyn.nc", "stsyn_points.nc", "stsyn_faults.nc", "stsyn_plain.nc"]:
    Path(name).unlink()

subprocess.run(["bash", "-ec", "source ../common.sh; remove_test_files"], check=True)
