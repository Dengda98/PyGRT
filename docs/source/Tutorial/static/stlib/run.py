# 有限断层输入文件由 run.sh 构建
import numpy as np
import pygrt

# BEGIN GRN
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="mod3")
pymod.static_greenfn(
    depsrc=np.arange(2, 6+1e-8, 1), deprcv=0.0, dists=np.arange(0, 36+1e-8, 2), calc_upar=True
)
# END GRN

# BEGIN SRC SYN
pymod.static_syn(
    norths=(-20, 20, 1), easts=(-20, 20, 1), deprcv=0.0,
    src_fault="fault.inp", src_fault_size=[1, 1], output_path="stsyn_ff.nc", zne=True, calc_upar=True
)
# END SRC SYN

# BEGIN COULOMB
pygrt.utils.static_stress("stsyn_ff.nc")
pygrt.utils.static_sproj("stsyn_ff.nc", strike=90, dip=60, rake=0)
pygrt.utils.static_coulomb("stsyn_ff.nc", friction=0.4)
# END COULOMB
