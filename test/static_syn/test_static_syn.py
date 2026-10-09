import subprocess
from pathlib import Path
from unittest import TestCase

import pygrt

subprocess.run(["bash", "-ec", "source ../common.sh; create_test_files"], check=True)

raises = TestCase().assertRaises

model = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="../milrow")
model.static_greenfn(depsrc=2, deprcv=0, dists=[0, 3, 6, 9, 12], calc_upar=True)
model.static_syn(scale=1e20, output_path="stsyn.nc")
model.static_syn(scale=1e20, force=[2, -1, 4], output_path="stsyn.nc")
model.static_syn(scale=1e20, strike=33, dip=44, rake=55, output_path="stsyn.nc")
model.static_syn(scale=1e6, scale_with_mu=True, strike=33, dip=44, output_path="stsyn.nc")
model.static_syn(scale=1e20, moment_tensor=[1, -2, -5, 0.5, 3, 1.2], zne=True, calc_upar=True, output_path="stsyn.nc")
model.static_syn(scale=1e20, norths=[-2, 2, 2], easts=[-2, 2, 2], output_path="stsyn.nc")

model = pygrt.PyModel1D(stgrn="stgrn_multi.nc", modelpath="../milrow")
model.static_greenfn(depsrc=[1, 2], deprcv=[0, 1], dists=[0, 3, 6, 9, 12], calc_upar=True)
model.static_syn(scale=1e20, depsrc=1.5, deprcv=0.5, norths=[-2, 2, 2], easts=[-2, 2, 2], nthreads=2, output_path="stsyn.nc")
model.static_syn(scale=1e20, depsrc=2, rcv_points="rcv_points.txt", output_path="stsyn.nc")
model.static_syn(scale=1e20, depsrc=2, rcv_points="rcv_geometry.txt", zne=True, calc_upar=True, output_path="stsyn.nc")
model.static_syn(scale=1e20, depsrc=2, rcv_fault="rcv_faults.inr", output_path="stsyn.nc")
model.static_syn(scale=1e20, depsrc=2, rcv_fault="rcv_faults.inr", rcv_fault_size=[1, 1], output_path="stsyn.nc")
model.static_syn(src_fault="faults.inp", src_fault_size=[1, 1], deprcv=0, calc_upar=True, output_path="stsyn.nc")
model.static_syn(src_fault="faults.inr", rcv_points="rcv_points.txt", output_path="stsyn.nc")

# 未指定点源源强
with raises(ValueError):
    model.static_syn(depsrc=2, deprcv=0, output_path="stsyn.nc")
# 使用接收点文件时，不可再设置接收深度
with raises(ValueError):
    model.static_syn(scale=1e20, depsrc=2, deprcv=0, rcv_points="rcv_points.txt", output_path="stsyn.nc")
# 接收断层的剖分尺寸不可为零或负数
with raises(ValueError):
    model.static_syn(scale=1e20, depsrc=2, rcv_fault="rcv_faults.inr", rcv_fault_size=[0, 1], output_path="stsyn.nc")

# 未指定多深度库对应的震源深度和接收深度
with raises(RuntimeError):
    model.static_syn(scale=1e20, output_path="stsyn.nc")
# 震源深度不可超出格林函数库范围
with raises(RuntimeError):
    model.static_syn(scale=1e20, depsrc=3, deprcv=0, output_path="stsyn.nc")
# 不可同时设置有限震源文件和点源源强
with raises(ValueError):
    model.static_syn(scale=1e20, src_fault="faults.inp", output_path="stsyn.nc")
# 不可同时设置接收点文件和接收网格
with raises(ValueError):
    model.static_syn(scale=1e20, depsrc=2, rcv_points="rcv_points.txt", norths=[-2, 2, 2], easts=[-2, 2, 2], output_path="stsyn.nc")

for name in ["stgrn.nc", "stgrn_multi.nc", "stsyn.nc"]:
    Path(name).unlink()

subprocess.run(["bash", "-ec", "source ../common.sh; remove_test_files"], check=True)
