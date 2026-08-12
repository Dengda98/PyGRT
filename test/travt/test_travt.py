import pygrt

depsrc = 2.0
deprcv = 0.0
modname = "../milrow"

pymod = pygrt.PyModel1D(modname)

for dist in [2, 3, 4, 5]:
    tp, ts = pymod.compute_travt1d(depsrc=depsrc, deprcv=deprcv, distarr=dist)
    print(dist, tp, ts)

# 多个震中距一次计算，分别返回 Tp、Ts 数组
tp, ts = pymod.compute_travt1d(depsrc=depsrc, deprcv=deprcv, distarr=[2, 3, 4, 5])
print(tp)
print(ts)
