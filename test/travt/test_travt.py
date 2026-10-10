from unittest import TestCase

import pygrt

raises = TestCase().assertRaises

model = pygrt.PyModel1D(modelpath="../milrow")
print(model.travt(depsrc=2, deprcv=0, dists=5))
print(model.travt(depsrc=2, deprcv=0, dists=[2, 5, 10]))

# 震中距列表不可包含逆序或重复值
with raises(ValueError):
    model.travt(depsrc=2, deprcv=0, dists=[5, 2])
# 震源深度不可为负数
with raises(ValueError):
    model.travt(depsrc=-1, deprcv=0, dists=5)
