from unittest import TestCase

import numpy as np
import pygrt

raises = TestCase().assertRaises

tbar = np.linspace(0, 2, 41)
pygrt.utils.lamb2(nu=0.25, tbar=tbar, R=10, depsrc=5, azimuth=30)
pygrt.utils.lamb2(nu=0.25, tbar=tbar, R=10, deprcv=5, azimuth=30)
pygrt.utils.lamb2(nu=0.25, tbar=tbar, R=10, depsrc=5, azimuth=30, phases=["P", "S", "SP"])
pygrt.utils.lamb2(nu=0.25, tbar=tbar, R=10, deprcv=5, azimuth=30, phases=["P", "S", "PS"])

# 第二类 Lamb 问题不可同时设置地下源和地下接收点
with raises(ValueError):
    pygrt.utils.lamb2(nu=0.25, tbar=tbar, R=10, depsrc=5, deprcv=3, azimuth=30)
# 地下源的深度不可为零或负数
with raises(ValueError):
    pygrt.utils.lamb2(nu=0.25, tbar=tbar, R=10, depsrc=0, azimuth=30)
# 水平距离不可为零或负数
with raises(ValueError):
    pygrt.utils.lamb2(nu=0.25, tbar=tbar, R=-1, depsrc=5, azimuth=30)
