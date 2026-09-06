import pygrt

# BEGIN LIBRARY
pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")
# eigenv 使用等间隔频率，供 modsum 进行逆傅里叶变换
pymod.eigenv(wtype="R", freqs=(0.0, 1.0, 0.02), phase_path="phase_R.nc", all_modes=True)
pymod.eigenv(wtype="L", freqs=(0.0, 1.0, 0.02), phase_path="phase_L.nc", all_modes=True)

# depsrc/deprcv 和 dists 传入列表，modsum 在内部遍历全部组合
# 加上 calc_upar=True 表示计算位移格林函数的空间偏导
for wave in ("R", "L"):
    pymod.modsum(
        phase_path=f"phase_{wave}.nc",
        depsrc=[2.0, 4.0],
        deprcv=[0.0, 2.0],
        dists=[80.0, 100.0, 120.0],
        modes=0,
        upsampling_n=4,
        calc_upar=True,
    )
# END LIBRARY


# BEGIN SYN
pymod = pygrt.PyModel1D(grn="GRN")
# syn 从面波格林函数库根目录精确选择一个深度和距离
# 加上 calc_upar=True 表示合成位移的空间偏导
pymod.syn(
    depsrc=4.0, deprcv=2.0, dist=100.0,
    azimuth=30.0, scale=1e24,
    strike=33.0, dip=50.0, rake=120.0,
    output_path="syn", calc_upar=True,
)
# END SYN
