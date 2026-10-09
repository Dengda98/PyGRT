from unittest import TestCase

import numpy as np
import pygrt

raises = TestCase().assertRaises

tbar = np.linspace(0, 2, 41)
pygrt.utils.lamb1(nu=0.25, tbar=tbar, azimuth=30)
pygrt.utils.lamb1(nu=0.25, tbar=tbar, azimuth=30, phases=["P", "S"])
pygrt.utils.lamb1(nu=0.25, tbar=tbar, azimuth=30, phases="R")
pygrt.utils.lamb1(nu=0.25, tbar=tbar, azimuth=30, cbar=2e-4)

# 泊松比不可超出 (0, 0.5) 范围
with raises(ValueError):
    pygrt.utils.lamb1(nu=0.5, tbar=tbar, azimuth=30)
# 移动源的接收点不可位于 x1 轴上
with raises(ValueError):
    pygrt.utils.lamb1(nu=0.25, tbar=tbar, azimuth=0, cbar=2e-4)

# 移动源的无量纲速度不可为零或负数
with raises(ValueError):
    pygrt.utils.lamb1(nu=0.25, tbar=tbar, azimuth=30, cbar=0)
# 相位列表不可为空
with raises(ValueError):
    pygrt.utils.lamb1(nu=0.25, tbar=tbar, azimuth=30, phases=[])
