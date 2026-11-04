import shutil
from pathlib import Path
from unittest import TestCase

import pygrt
from pygrt import noise

raises = TestCase().assertRaises
model = pygrt.PyModel1D(modelpath="../milrow")
model.synnoise(source_file="noise_sources.txt", receiver_file="noise_receivers.txt",
               output_path="noise_python", modes=None, nthreads=2, print_log=False)

points = noise.random_points(2, radius_range=(2, 3), depth=0.5, seed=0)
noise.write_sources("helper_sources.txt", points, interval=0.5, force_z_range=(-0.5, 1),
                    force_azimuth_range=(30, 150), strength_range=(1e15, 3e15), seed=0)
noise.write_receivers("helper_receivers.txt", [[0, 0, 0]], dt=0.2, length=2, ids=["A"])
model.synnoise(source_file="helper_sources.txt", receiver_file="helper_receivers.txt",
               output_path="noise_helper", nthreads=1, print_log=False)
noise.random_points(2, north_range=(-1, 1), east_range=(2, 3), seed=0)

# 必须指定介质模型
with raises(RuntimeError):
    pygrt.PyModel1D().synnoise(source_file="noise_sources.txt", receiver_file="noise_receivers.txt",
                              output_path="invalid")
# 只支持自由表面和半空间底边界
with raises(ValueError):
    pygrt.PyModel1D(modelpath="../milrow", topbound="rigid").synnoise(
        source_file="noise_sources.txt", receiver_file="noise_receivers.txt", output_path="invalid")
# 不可混用矩形与环带范围
with raises(ValueError):
    noise.random_points(2, radius_range=(2, 3), north_range=(-1, 1))
# 激发间隔必须为正
with raises(ValueError):
    noise.write_sources("invalid_sources.txt", points, interval=0)
# 采样间隔必须为正
with raises(ValueError):
    noise.write_receivers("invalid_receivers.txt", [[0, 0, 0]], dt=0, length=2)
# 记录长度必须为正
with raises(ValueError):
    noise.write_receivers("invalid_receivers.txt", [[0, 0, 0]], dt=0.2, length=0)
# 记录长度必须对应整数采样点
with raises(ValueError):
    noise.write_receivers("invalid_receivers.txt", [[0, 0, 0]], dt=0.2, length=2.1)

shutil.rmtree("noise_python")
shutil.rmtree("noise_helper")
Path("helper_sources.txt").unlink()
Path("helper_receivers.txt").unlink()
