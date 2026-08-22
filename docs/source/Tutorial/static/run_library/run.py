import pygrt

# BEGIN GRN
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")
# 加上 calc_upar=True 表示计算位移格林函数的空间偏导
pymod.compute_static_grn(
    depsrc=[2.0, 4.0], 
    deprcv=[0.0, 2.0], 
    dists=[0.0, 5.0, 10.0, 15.0], 
    calc_upar=True
)
# END GRN


# BEGIN SYN
# 指定静态格林函数路径
pymod = pygrt.PyModel1D(stgrn="stgrn.nc")

# 加上 calc_upar=True 表示合成位移的空间偏导
# 需要设置 depsrc, deprcv 以明确指出震源深度和台站深度
pymod.compute_static_syn(
    norths=(-10.0, 10.0, 5.0), easts=(-10.0, 10.0, 5.0),
    depsrc=3.0, deprcv=1.0,
    scale=1e24, strike=33.0, dip=50.0, rake=120.0, output_path="stsyn.nc", calc_upar=True
)
# END SYN
