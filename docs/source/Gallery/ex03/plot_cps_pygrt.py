import pygrt
import numpy as np 
import matplotlib.pyplot as plt 
from obspy import * 


depsrc = 2.0   
deprcv = 0     
pymod = pygrt.PyModel1D(grn="GRN_pygrt", modelpath="./milrow")

rs = np.array([10]) 

nt = 1024   
dt = 0.01   
zeta = 0.8  

# compute Green's Functions
pymod.compute_grn(
    depsrc=depsrc,
    deprcv=deprcv,
    dists=rs,
    nt=nt,
    dt=dt,
    zeta=zeta,
    Length=20,
)
# 仅一个震中距，可用通配符读回
st_grt = read("GRN_pygrt/*/*.sac")


try:
    st_cps = read("milrow_sdep2_rdep0/GRN/*")

    from plot_cps_grt import plot
    plot(st_grt, st_cps, "compare_cps_pygrt.svg")
except Exception as e:
    print(str(e))

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
p = Path("GRN_pygrt")
if p.is_dir():
    shutil.rmtree(p, ignore_errors=True)
