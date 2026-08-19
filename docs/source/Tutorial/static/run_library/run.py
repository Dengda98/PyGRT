from pathlib import Path

import pygrt


MODEL = Path("milrow")


# BEGIN GRN
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=MODEL)
# 加上 calc_upar=True 表示计算位移格林函数的空间偏导
pymod.compute_static_grn(depsrc=[2.0, 4.0], deprcv=[0.0, 2.0], dists=[0.0, 5.0, 10.0, 15.0], calc_upar=True)
# END GRN


# BEGIN SYN
pymod = pygrt.PyModel1D(stgrn="stgrn.nc")
# 加上 calc_upar=True 表示合成位移的空间偏导
pymod.compute_static_syn(
    scale=1e24,
    output_path="stsyn.nc",
    strike=33.0,
    dip=50.0,
    rake=120.0,
    depsrc=3.0,
    deprcv=1.0,
    norths=(-10.0, 10.0, 5.0),
    easts=(-10.0, 10.0, 5.0),
    calc_upar=True
)
# END SYN
