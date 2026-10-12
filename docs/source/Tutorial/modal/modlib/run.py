import pygrt


# BEGIN LIBRARY
pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")
# eigenv 使用等间隔频率，供 modsum 进行逆傅里叶变换
pymod.eigenv(wtype="R", freqs=(0.0, 1.0, 0.02), phase_path="phase_R.nc", all_modes=True)
pymod.eigenv(wtype="L", freqs=(0.0, 1.0, 0.02), phase_path="phase_L.nc", all_modes=True)

# 计算全部震源深度、台站深度和震中距的组合
for wave in ("R", "L"):
    pymod.modsum(
        phase_path=f"phase_{wave}.nc",
        depsrc=[2.0, 4.0], deprcv=[0.0, 2.0], dists=[80.0, 100.0, 120.0],
        modes=0, upsampling_n=4,
    )
# END LIBRARY
