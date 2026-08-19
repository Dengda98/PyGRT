from pathlib import Path

import pygrt


MODEL = Path("milrow")


# BEGIN GRN
pymod = pygrt.PyModel1D(grn="GRN", modelpath=MODEL)
# 加上 calc_upar=True 表示计算位移格林函数的空间偏导
pymod.compute_grn(depsrc=[2.0, 4.0], deprcv=[0.0, 2.0], dists=[5.0, 8.0, 10.0], nt=256, dt=0.02, calc_upar=True)
# END GRN


# BEGIN SYN
pymod = pygrt.PyModel1D(grn="GRN")
# 加上 calc_upar=True 表示合成位移的空间偏导
pymod.compute_syn(depsrc=4.0, deprcv=2.0, dist=8.0, azimuth=30.0, scale=1e24, output_path="syn_python", calc_upar=True)
# END SYN
