import shutil
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

medium = [8, 4.62, 3.3]
dt = 0.1
pygrt.utils.lamb(modelparams=medium, nt=32, dt=dt, dist=10, depsrc=0, deprcv=0,
                 azimuth=30, force=[0.5, -1, 2], scale=1e15, output_path="lamb")
pygrt.utils.lamb(modelparams=medium, nt=32, dt=dt, dist=10, depsrc=5, deprcv=0,
                 azimuth=30, strike=100, dip=30, rake=70, scale=1e20, calc_upar=True, zne=True, output_path="lamb")
pygrt.utils.lamb(modelparams=medium, nt=32, dt=dt, dist=10, depsrc=0, deprcv=5,
                 azimuth=30, moment_tensor=[1, 0, 0, 1, 0, 1], scale=1e20, output_path="lamb")
pygrt.utils.lamb(modelparams=medium, nt=32, dt=dt, dist=10, depsrc=5, deprcv=1,
                 azimuth=30, strike=100, dip=30, scale=1e10, scale_with_mu=True, output_path="lamb")
pygrt.utils.lamb(modelparams=medium, nt=32, dt=dt, dist=10, depsrc=5, deprcv=0,
                 azimuth=30, scale=1e20, phases=["P", "S", "SP"], output_path="lamb")
pygrt.utils.lamb(modelparams=medium, nt=32, dt=dt, dist=10, depsrc=5, deprcv=1, azimuth=30, scale=1e20,
                 time_function="p/0.5", integrate_order=1, differentiate_order=1, delayT0=0.2, delayV0=9, output_path="lamb")
pygrt.utils.lamb(modelparams=medium, nt=32, dt=dt, dist=10, depsrc=5, deprcv=1,
                 azimuth=30, scale=1e20, ref_first_p=True, delayT0=-1, output_path="lamb")
# 自定义时间函数按 dt 采样
pygrt.utils.lamb(modelparams=medium, nt=32, dt=dt, dist=10, depsrc=5, deprcv=1,
                 azimuth=30, scale=1e20, time_function="0/time_function.txt+d0.2", output_path="lamb")
# 两列时间、振幅输入支持非等距采样
pygrt.utils.lamb(modelparams=medium, nt=32, dt=dt, dist=10, depsrc=5, deprcv=1,
                 azimuth=30, scale=1e20, time_function="0/time_function_nonuniform.txt+d0.2", output_path="lamb")
pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, depsrc=5, scale=1e20,
                 rcv_points="rcv_geometry.txt", calc_upar=True, nthreads=2, output_path="lamb_points")
pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, depsrc=5, scale=1e20,
                 rcv_fault="rcv_faults.inr", rcv_fault_size=[1, 1], output_path="lamb_faults")
pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, src_fault="faults.inp", src_fault_size=[1, 1],
                 rcv_points="rcv_points.txt", print_log=False, output_path="lamb_points")
pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, src_fault="faults.inr", src_fault_size=[1, 1],
                 rcv_fault="rcv_faults.inr", print_log=False, output_path="lamb_faults")

# 显式不再剖分，每条矩形源及接收断层只取一个中心点
pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, src_fault="faults.inp", src_fault_size=(0, 0),
                 rcv_fault="rcv_faults.inr", rcv_fault_size=(0, 0), print_log=False, output_path="lamb_faults")

# 各断层的破裂过程和延迟，以及全局时间函数覆盖行末设置
pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, src_fault="faults_rupture.inr", src_fault_size=[1, 1],
                 rcv_points="rcv_points.txt", print_log=False, output_path="lamb_points")
pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, src_fault="faults_rupture.inr", src_fault_size=[1, 1],
                 rcv_points="rcv_points.txt", time_function="t/0.2/0.1/0.3+d0.1", print_log=False, output_path="lamb_points")

# 有限断层的时间函数不可只设置在部分记录中
with raises(RuntimeError):
    pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, src_fault="faults_partial.inr", src_fault_size=[1, 1],
                     rcv_points="rcv_points.txt", output_path="lamb_points")
# 破裂延迟不可为负数
with raises(RuntimeError):
    pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, src_fault="faults_negative_delay.inr", src_fault_size=[1, 1],
                     rcv_points="rcv_points.txt", output_path="lamb_points")

# 剖分尺寸不可为负数
with raises(ValueError):
    pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, src_fault="faults.inp", src_fault_size=(-1, -1),
                     rcv_points="rcv_points.txt", output_path="lamb_points")

# 未指定点源源强
with raises(ValueError):
    pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, dist=10, depsrc=5, deprcv=1, azimuth=30, output_path="lamb")
# 矩形有限震源不可省略剖分尺寸
with raises(RuntimeError):
    pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, src_fault="faults.inr", rcv_points="rcv_points.txt", output_path="lamb_points")

# 源和接收点都在地表时，不可使用爆炸源
with raises(RuntimeError):
    pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, dist=10, depsrc=0, deprcv=0, azimuth=30, scale=1e20, output_path="lamb")
# 水平距离不可为零或负数
with raises(RuntimeError):
    pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, dist=0, depsrc=5, deprcv=1, azimuth=30, scale=1e20, output_path="lamb")
# 不可同时设置接收点文件和接收断层文件
with raises(ValueError):
    pygrt.utils.lamb(modelparams=medium, nt=16, dt=dt, depsrc=5, scale=1e20,
                     rcv_points="rcv_points.txt", rcv_fault="rcv_faults.inr", output_path="lamb")

for name in ["lamb", "lamb_points", "lamb_faults"]:
    shutil.rmtree(name)
