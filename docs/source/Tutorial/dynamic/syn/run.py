import pygrt

# 独立运行本页时，先计算与上一节相同的格林函数
pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")
pymod.greenfn(depsrc=2.0, deprcv=0.0, dists=[5, 8, 10], nt=500, dt=0.02)

# BEGIN REUSE GRN
# 若仅使用已算好的格林函数做合成，构造时只需指定 grn，无需 modelpath
pymod = pygrt.PyModel1D(grn="GRN")
# END REUSE GRN

# BEGIN SYN SOURCES
# 爆炸源，标量矩 1e24 dyne·cm
pymod.syn(dist=10.0, azimuth=30.0, scale=1e24, output_path="syn_ex")

# 单力源，(fN, fE, fZ)=(1, -0.5, 2)，标度 1e16 dyne
pymod.syn(
    dist=10.0, azimuth=30.0, scale=1e16, output_path="syn_sf",
    force=(1, -0.5, 2),
)

# 剪切源，strike=33°，dip=50°，rake=120°，标量矩 1e24 dyne·cm
pymod.syn(
    dist=10.0, azimuth=30.0, scale=1e24, output_path="syn_dc",
    strike=33, dip=50, rake=120,
)

# 走滑剪切源，strike=33°，dip=90°，rake=0°，标量矩 1e24 dyne·cm
pymod.syn(
    dist=10.0, azimuth=30.0, scale=1e24, output_path="syn_dc2",
    strike=33, dip=90, rake=0,
)

# 张裂源，strike=33°，dip=50°，标量矩 1e24 dyne·cm
pymod.syn(
    dist=10.0, azimuth=30.0, scale=1e24, output_path="syn_ts",
    strike=33, dip=50,
)

# 矩张量源，标量矩 1e24 dyne·cm
pymod.syn(
    dist=10.0, azimuth=30.0, scale=1e24, output_path="syn_mt",
    moment_tensor=(0.1, -0.2, 1.0, 0.3, -0.5, -2.0),
)
# END SYN SOURCES

# BEGIN ZNE
# 接之前的代码
# 设置 zne=True 可返回 ZNE 分量
pymod.syn(dist=10.0, azimuth=30.0, scale=1e24, output_path="syn_dc_zne", strike=33, dip=50, rake=120, zne=True)
# END ZNE

# BEGIN TIME FUNC
# time_function 对应 CLI -D
# 本例中的震源时间函数均使用面积归一化
# 平台段时长为零时梯形波退化为三角波
# 卷积用的时间函数会以 sig.sac 保存在输出目录
pymod.syn(dist=10.0, azimuth=30.0, scale=1e16, output_path="syn_sf_trig", force=(1, -0.5, 2), time_function="t/0.3/0/0.3")
# 非对称余弦波，上升段时长/下降段时长
pymod.syn(dist=10.0, azimuth=30.0, scale=1e16, output_path="syn_sf_cosine", force=(1, -0.5, 2), time_function="c/0.1/0.5")
# END TIME FUNC

# BEGIN INT DIF
# 在合成时指定积分或微分次数，分别对应 CLI 的 -I 和 -J
pymod.syn(dist=10.0, azimuth=30.0, scale=1e24, output_path="syn_mt_int1", moment_tensor=(0.1, -0.2, 1.0, 0.3, -0.5, -2.0), integrate_order=1)
pymod.syn(dist=10.0, azimuth=30.0, scale=1e24, output_path="syn_mt_dif1", moment_tensor=(0.1, -0.2, 1.0, 0.3, -0.5, -2.0), differentiate_order=1)
# END INT DIF
