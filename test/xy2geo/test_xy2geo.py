from pathlib import Path
from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

pygrt.utils.xy2geo(qfile="rcv_points.txt", outgrid="geo.txt", lat0=35, lon0=179.5)
pygrt.utils.geo2xy(qfile="geo.txt", outgrid="local.txt", lat0=35, lon0=179.5)

medium = [6, 3.464, 2.7]
pygrt.utils.okada(modelparams=medium, scale=1e20, depsrc=10, deprcv=0,
                  norths=[-2, 2, 2], easts=[-2, 2, 2], output_path="grid.nc")
pygrt.utils.xy2geo("grid.nc", outgrid="geo.nc", lat0=35, lon0=179.5)
pygrt.utils.geo2xy("geo.nc", outgrid="local.nc", lat0=35, lon0=179.5)
pygrt.utils.okada(modelparams=medium, scale=1e20, depsrc=10, rcv_points="rcv_points.txt", output_path="points.nc")
pygrt.utils.xy2geo("points.nc", outgrid="geo.nc", lat0=-20, lon0=-179.5)
pygrt.utils.geo2xy("geo.nc", outgrid="local.nc", lat0=-20, lon0=-179.5)
pygrt.utils.okada(modelparams=medium, scale=1e20, depsrc=10, rcv_fault="rcv_faults.inr", output_path="faults.nc")
pygrt.utils.xy2geo("faults.nc", outgrid="geo.nc", lat0=40, lon0=10)
pygrt.utils.geo2xy("geo.nc", outgrid="local.nc", lat0=40, lon0=10)

# 输入和输出不可使用同一个文件
with raises(RuntimeError):
    pygrt.utils.xy2geo(qfile="rcv_points.txt", outgrid="rcv_points.txt", lat0=35, lon0=10)

# 参考经度不可超出 [-180, 180] 度
with raises(ValueError):
    pygrt.utils.xy2geo(qfile="rcv_points.txt", outgrid="geo.txt", lat0=35, lon0=181)

# 不可同时设置 NetCDF 文件和文本坐标文件
with raises(ValueError):
    pygrt.utils.xy2geo("grid.nc", qfile="rcv_points.txt", outgrid="geo.nc", lat0=35, lon0=10)
# 参考纬度不可达到或超过 ±90 度
with raises(ValueError):
    pygrt.utils.xy2geo(qfile="rcv_points.txt", outgrid="geo.txt", lat0=91, lon0=10)
# 未指定参考经纬度
with raises(ValueError):
    pygrt.utils.geo2xy(qfile="geo.txt", outgrid="local.txt")

for name in ["grid.nc", "points.nc", "faults.nc", "geo.nc", "local.nc", "geo.txt", "local.txt"]:
    Path(name).unlink()
