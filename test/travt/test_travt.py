import pygrt

depsrc = 2.0
deprcv = 0.0
modname = "../milrow"

pymod = pygrt.PyModel1D(modelpath=modname)

for dist in [2, 3, 4, 5]:
    tp, ts = pymod.compute_travt1d(depsrc=depsrc, deprcv=deprcv, dists=dist)
    print(dist, tp, ts)

# 多个震中距一次计算，分别返回 Tp、Ts 数组
tp, ts = pymod.compute_travt1d(depsrc=depsrc, deprcv=deprcv, dists=[2, 3, 4, 5])
print(tp)
print(ts)

try:
    pymod.compute_travt1d(depsrc=depsrc, deprcv=deprcv, dists=[2, 1])
except ValueError:
    pass
else:
    raise AssertionError("non-ascending dists should raise")
