# 有限断层输入文件由 run.sh 构建
import shutil

import numpy as np
import pygrt

# 重新建库，避免保留旧的接收深度采样节点
shutil.rmtree("GRN", ignore_errors=True)

# BEGIN GRN
pymod = pygrt.PyModel1D(grn="GRN", modelpath="mod3")
pymod.greenfn(
    depsrc=np.arange(2, 6+1e-8, 1), deprcv=0.0, dists=np.arange(2, 36+1e-8, 2),
    nt=512, dt=0.1, keepAllFreq=True, safilonTol=1e-3, calc_upar=True
)
# END GRN

# BEGIN SRC SYN
# 每行已是子断层，指定不再剖分以保留每行一个子源；不覆盖行末时间函数
pymod.syn(src_fault="fault_rupture.inp", src_fault_size=(0, 0), dist=20.0, azimuth=30.0, deprcv=0.0,
          output_path="syn_ff", calc_upar=True, integrate_order=1)
pygrt.utils.stress("syn_ff")
pygrt.utils.sproj("syn_ff", strike=90, dip=60, rake=0)
pygrt.utils.coulomb("syn_ff", friction=0.4)
# END SRC SYN
