from pathlib import Path
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

medium = [6, 3.464, 2.7]
pygrt.utils.okada(modelparams=medium, scale=1e12, scale_with_mu=True, depsrc=10, deprcv=0,
                  norths=[-2, 2, 2], easts=[-2, 2, 2], output_path="okada.nc")
pygrt.utils.okada(modelparams=medium, scale=1e16, scale_with_mu=True, depsrc=10, deprcv=0,
                  norths=[-2, 2, 2], easts=[-2, 2, 2], strike=100, dip=20, rake=80, zne=True, calc_upar=True, output_path="okada.nc")
pygrt.utils.okada(modelparams=medium, scale=1e16, depsrc=10, deprcv=0,
                  norths=[-2, 2, 2], easts=[-2, 2, 2], strike=100, dip=20, output_path="okada.nc")
pygrt.utils.okada(modelparams=medium, scale=1e12, depsrc=10, rcv_points="rcv_geometry.txt", output_path="okada.nc")
pygrt.utils.okada(modelparams=medium, scale=1e12, depsrc=10, rcv_fault="rcv_faults.inr", output_path="okada.nc")
pygrt.utils.okada(modelparams=medium, scale=1e12, depsrc=10, rcv_fault="rcv_faults.inr", rcv_fault_size=[1, 1], output_path="okada.nc")
pygrt.utils.okada(modelparams=medium, src_fault="faults.inp", deprcv=0,
                  norths=[-2, 2, 2], easts=[-2, 2, 2], zne=True, calc_upar=True, output_path="okada.nc")
pygrt.utils.okada(modelparams=medium, src_fault="faults.inr", rcv_points="rcv_points.txt", output_path="okada.nc")

# 接收网格不可只设置北向坐标
with raises(ValueError):
    pygrt.utils.okada(modelparams=medium, scale=1e20, depsrc=10, deprcv=0, norths=[-2, 2, 2], output_path="okada.nc")
# 使用接收点文件时，不可再设置接收深度
with raises(ValueError):
    pygrt.utils.okada(modelparams=medium, scale=1e20, depsrc=10, deprcv=0, rcv_points="rcv_points.txt", output_path="okada.nc")

# 震源深度不可为负数
with raises(ValueError):
    pygrt.utils.okada(modelparams=medium, scale=1e12, depsrc=-1, rcv_points="rcv_points.txt", output_path="okada.nc")
# 不可同时设置有限震源文件和点源源强
with raises(ValueError):
    pygrt.utils.okada(modelparams=medium, scale=1e12, src_fault="faults.inp", rcv_points="rcv_points.txt", output_path="okada.nc")
# 不可同时设置接收点文件和接收断层文件
with raises(ValueError):
    pygrt.utils.okada(modelparams=medium, scale=1e12, depsrc=10, rcv_points="rcv_points.txt",
                      rcv_fault="rcv_faults.inr", output_path="okada.nc")

Path("okada.nc").unlink()
