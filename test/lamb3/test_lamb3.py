from unittest import TestCase

import numpy as np
import pygrt

raises = TestCase().assertRaises

tbar = np.linspace(0, 2, 41)
pygrt.utils.lamb3(nu=0.25, tbar=tbar, R=10, depsrc=2, deprcv=1, azimuth=30)
pygrt.utils.lamb3(nu=0.25, tbar=tbar, R=10, depsrc=2, deprcv=1, azimuth=30, phases=["P", "S", "PP", "SS", "PS", "SP", "sPs"])

# 地下源的深度不可为零或负数
with raises(ValueError):
    pygrt.utils.lamb3(nu=0.25, tbar=tbar, R=10, depsrc=0, deprcv=1, azimuth=30)
# 地下接收点的深度不可为零或负数
with raises(ValueError):
    pygrt.utils.lamb3(nu=0.25, tbar=tbar, R=10, depsrc=2, deprcv=0, azimuth=30)
# 水平距离不可为零或负数
with raises(ValueError):
    pygrt.utils.lamb3(nu=0.25, tbar=tbar, R=0, depsrc=2, deprcv=1, azimuth=30)
