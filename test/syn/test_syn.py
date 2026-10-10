import shutil
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

dt = 0.1

model = pygrt.PyModel1D(grn="GRN", modelpath="../milrow")
model.greenfn(depsrc=2, deprcv=0, dists=5, nt=32, dt=dt, calc_upar=True, print_log=False)
model.syn(azimuth=22, scale=1e20, output_path="syn")
model.syn(azimuth=22, scale=1e16, force=[-1, 2, -4], output_path="syn")
model.syn(azimuth=22, scale=1e20, strike=33, dip=44, rake=55, output_path="syn")
model.syn(azimuth=22, scale=1e10, scale_with_mu=True, strike=33, dip=44, output_path="syn")
model.syn(azimuth=22, scale=1e20, moment_tensor=[1, -2, -5, 0.5, 3, 1.2], output_path="syn")
model.syn(azimuth=22, scale=1e20, time_function="p/0.6", output_path="syn")
model.syn(azimuth=22, scale=1e20, time_function="t/0.2/0.2/0.3", output_path="syn")
model.syn(azimuth=22, scale=1e20, time_function="c/0.2/0.4", output_path="syn")
model.syn(azimuth=22, scale=1e20, time_function="r/1.2", output_path="syn")
# 时间函数按格林函数的 dt 采样
model.syn(azimuth=22, scale=1e20, time_function="0/time_function.txt", output_path="syn")
# 两列输入从 0.0 开始，支持不同采样间隔和非等距采样
model.syn(azimuth=22, scale=1e20, time_function="0/time_function_times.txt", output_path="syn")
model.syn(azimuth=22, scale=1e20, time_function="0/time_function_nonuniform.txt+d0.2", output_path="syn")
model.syn(azimuth=22, scale=1e20, integrate_order=1, output_path="syn")
model.syn(azimuth=22, scale=1e20, differentiate_order=1, zne=True, calc_upar=True, output_path="syn")

node = pygrt.PyModel1D(grn="GRN/milrow_2_0_5", modelpath="../milrow")
node.syn(azimuth=22, scale=1e20, output_path="syn")
model = pygrt.PyModel1D(grn="GRN_MULTI", modelpath="../milrow")
model.greenfn(depsrc=[1, 2], deprcv=[0, 1], dists=[1, 6, 12], nt=16, dt=dt, calc_upar=True, print_log=False)
model.syn(depsrc=1.5, deprcv=0.5, dist=8, azimuth=22, scale=1e20, output_path="syn")
model.syn(depsrc=1.5, deprcv=0.5, dist=8, azimuth=22, scale=1e20, interpolate=False, output_path="syn")
model.syn(depsrc=2, scale=1e20, rcv_points="rcv_geometry.txt", nthreads=2, calc_upar=True, output_path="syn_points")
model.syn(depsrc=2, scale=1e20, rcv_fault="rcv_faults.inr", output_path="syn_faults")
model.syn(src_fault="faults.inp", src_fault_size=[1, 1], rcv_points="rcv_points.txt", output_path="syn_points")
model.syn(src_fault="faults.inr", rcv_fault="rcv_faults.inr", rcv_fault_size=[1, 1], output_path="syn_faults")

# 显式不再剖分，每条矩形源及接收断层只取一个中心点
model.syn(src_fault="faults.inp", src_fault_size=(0, 0), rcv_fault="rcv_faults.inr",
          rcv_fault_size=(0, 0), calc_upar=True, output_path="syn_faults")

# 各断层的破裂过程和延迟，以及全局时间函数覆盖行末设置
model.syn(src_fault="faults_rupture.inr", src_fault_size=[1, 1], rcv_points="rcv_points.txt", output_path="syn_points")
model.syn(src_fault="faults_rupture.inr", src_fault_size=[1, 1], rcv_points="rcv_points.txt",
          time_function="p/0.4+d0.1", output_path="syn_points")
model.syn(src_fault="faults_rupture.inr", src_fault_size=[1, 1], rcv_points="rcv_points.txt",
          time_function="0/time_function_nonuniform.txt", output_path="syn_points")

# 有限断层的时间函数不可只设置在部分记录中
with raises(RuntimeError):
    model.syn(src_fault="faults_partial.inr", src_fault_size=[1, 1], rcv_points="rcv_points.txt", output_path="syn_points")
# 破裂延迟不可为负数
with raises(RuntimeError):
    model.syn(src_fault="faults_negative_delay.inr", src_fault_size=[1, 1], rcv_points="rcv_points.txt", output_path="syn_points")

# 剖分尺寸须同时为零或同时为正数
with raises(ValueError):
    model.syn(src_fault="faults.inp", src_fault_size=(0, 1), rcv_points="rcv_points.txt", output_path="syn_points")

# 未指定点源源强
with raises(ValueError):
    model.syn(depsrc=2, deprcv=0, dist=6, azimuth=22, output_path="syn")
# 使用接收点文件时，不可再设置接收深度
with raises(ValueError):
    model.syn(depsrc=2, deprcv=0, scale=1e20, rcv_points="rcv_points.txt", output_path="syn")
# 使用格林函数节点子目录时，不可再设置深度和距离选择参数
with raises(RuntimeError):
    node.syn(depsrc=2, azimuth=22, scale=1e20, output_path="syn")
# 抛物线时间函数的持续时间不可为零或负数
with raises(RuntimeError):
    node.syn(azimuth=22, scale=1e20, time_function="p/0", output_path="syn")

# 未指定多深度、多距离库对应的震源深度、接收深度和距离
with raises(RuntimeError):
    model.syn(azimuth=22, scale=1e20, output_path="syn")
# 震中距不可超出格林函数库范围
with raises(RuntimeError):
    model.syn(depsrc=2, deprcv=0, dist=13, azimuth=22, scale=1e20, output_path="syn")
# 不可同时设置单力源和双力偶源机制
with raises(ValueError):
    model.syn(depsrc=2, scale=1e20, force=[1, 2, 3], strike=33, dip=44, output_path="syn")
# 不可同时设置接收点文件和接收断层文件
with raises(ValueError):
    model.syn(depsrc=2, scale=1e20, rcv_points="rcv_points.txt", rcv_fault="rcv_faults.inr", output_path="syn")

for name in ["GRN", "GRN_MULTI", "syn", "syn_points", "syn_faults"]:
    shutil.rmtree(name)
