import numpy as np
import pygrt 

pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")

pymod.greenfn(
    depsrc=5.0,
    deprcv=0.0,
    dists=[2500],
    nt=2000,
    dt=1,
    Length=20,
    safilonTol=1e-2,  # 自适应采样精度
    filonCut=10,
    delayT0=100,
)

# 删除中间计算结果（成图由 plot.py 负责）
import shutil
from pathlib import Path
p = Path("GRN")
if p.is_dir():
    shutil.rmtree(p, ignore_errors=True)
