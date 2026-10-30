# BEGIN GRN
import numpy as np
import pygrt
from obspy import read

pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")

# 传入 calc_upar=True 计算空间导数
pymod.greenfn(depsrc=2.0, deprcv=0.0, dists=[30], nt=1500, dt=0.02, calc_upar=True)
stgrn = read("GRN/*/*.sac")
print(stgrn.__str__(extended=True))
# 45 Trace(s) in Stream:
# .SYN..EXZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:29.980000Z | 50.0 Hz, 1500 samples
# ...
# .SYN..zEXZ | ...
# .SYN..rEXZ | ...
# END GRN

# BEGIN SYN DC
# 传入 calc_upar=True 计算空间导数，integrate_order=1 得到阶跃源响应
# time_function 卷积上升 0.5 s、下降 1.5 s 的非对称余弦波
# ?.sac 为位移，[zrt]?.sac 为空间导数
pymod.syn(dist=30.0, azimuth=30.0, scale=1e24, output_path="syn_dc",
          strike=33, dip=50, rake=120, calc_upar=True, integrate_order=1, time_function="c/0.5/1.5")
stsyn = read("syn_dc/?.sac") + read("syn_dc/[zrt]?.sac")
print(stsyn)
# 12 Trace(s) in Stream:
# .SYN..Z  | ...
# .SYN..zZ | ...
# .SYN..rZ | ...
# .SYN..tZ | ...
# END SYN DC

# BEGIN ZNE
# 传入 zne=True 可返回 ZNE 分量，integrate_order=1 得到阶跃源响应
# ?.sac 为位移，[zne]?.sac 为空间导数
pymod.syn(dist=30.0, azimuth=30.0, scale=1e24, output_path="syn_dc_zne",
          strike=33, dip=50, rake=120, calc_upar=True, zne=True, integrate_order=1, time_function="c/0.5/1.5")
stsyn = read("syn_dc_zne/?.sac") + read("syn_dc_zne/[zne]?.sac")
print(stsyn)
# 12 Trace(s) in Stream:
# .SYN..Z  | ...
# .SYN..zZ | ...
# .SYN..nZ | ...
# .SYN..eZ | ...
# END ZNE

# BEGIN STRAIN
pygrt.tensor.strain("syn_dc_zne")
st_strain = read("syn_dc_zne/strain_*.sac")
print(st_strain)
# 6 Trace(s) in Stream:
# .SYN..ZZ | ...
# .SYN..ZN | ...
# .SYN..ZE | ...
# .SYN..NN | ...
# .SYN..NE | ...
# .SYN..EE | ...
# END STRAIN

# BEGIN ROTATION
pygrt.tensor.rotation("syn_dc_zne")
st_rotation = read("syn_dc_zne/rotation_*.sac")
print(st_rotation)
# 3 Trace(s) in Stream:
# .SYN..ZN | ...
# .SYN..ZE | ...
# .SYN..NE | ...
# END ROTATION

# BEGIN STRESS
pygrt.tensor.stress("syn_dc_zne")
st_stress = read("syn_dc_zne/stress_*.sac")
print(st_stress)
# 6 Trace(s) in Stream:
# .SYN..ZZ | ...
# .SYN..ZN | ...
# .SYN..ZE | ...
# .SYN..NN | ...
# .SYN..NE | ...
# .SYN..EE | ...
# END STRESS

# BEGIN SPROJ
# 接收断层形态与震源机制独立设置
pygrt.tensor.sproj("syn_dc_zne", strike=90, dip=60, rake=0)
# END SPROJ

# BEGIN COULOMB
pygrt.tensor.coulomb("syn_dc_zne", friction=0.4)
# END COULOMB
