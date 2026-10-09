import pygrt


# BEGIN GRN
pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")
pymod.greenfn(depsrc=[2.0, 4.0], deprcv=[0.0, 2.0], dists=[5.0, 8.0, 10.0], nt=256, dt=0.02)
# END GRN


# BEGIN SYN
pymod = pygrt.PyModel1D(grn="GRN")
# 目标位置位于库范围内，无需与采样节点重合
# interpolate=True 使用线性插值，也是默认方式
pymod.syn(depsrc=3.0, deprcv=1.0, dist=7.0, azimuth=30.0, scale=1e24, interpolate=True, output_path="syn_linear")

# 同一目标位置，改用 interpolate=False 选择最近邻
pymod.syn(depsrc=3.0, deprcv=1.0, dist=7.0, azimuth=30.0, scale=1e24, interpolate=False, output_path="syn_nearest")
# END SYN
