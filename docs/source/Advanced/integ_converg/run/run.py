# -------------------------------------------------------------------
# BEGIN DGRN
import numpy as np
import pygrt 
from pygrt.cli import format_float

depsrc = 2.0
deprcv = 0.0

pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")

# statsidxs 指定频率索引，核函数写入 GRN_grtstats/{model}_{depsrc}_{deprcv}/
dists = [5,8,10]
pymod.compute_grn(depsrc=depsrc, deprcv=deprcv, dists=dists, nt=500, dt=0.02, statsidxs=[50,100])
# END DGRN
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# BEGIN read statsfile
# 可使用通配符简化输入，因为对应索引值下只会有一个文件
# 返回的是自定义类型的numpy数组
statsdir = f"GRN_grtstats/milrow_{format_float(depsrc)}_{format_float(deprcv)}"
statsdata = pygrt.utils.read_statsfile(f"{statsdir}/K_0050_*")
print(statsdata.dtype)
# [('k', '<f8'), ('EX_q', '<c16'), ('EX_w', '<c16'), ('VF_q', '<c16'), ('VF_w', '<c16'), ('HF_q', '<c16'), ('HF_w', '<c16'), ('HF_v', '<c16'), ('DD_q', '<c16'), ('DD_w', '<c16'), ('DS_q', '<c16'), ('DS_w', '<c16'), ('DS_v', '<c16'), ('SS_q', '<c16'), ('SS_w', '<c16'), ('SS_v', '<c16')]
# END read statsfile
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# BEGIN plot stats
ir = 2
dist=dists[ir]
srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata(statsdata, dist=dist, srctype=srctype, ptype=ptype)
fig.savefig(f"{srctype}_{ptype}.svg", bbox_inches='tight')
# END plot stats
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# BEGIN plot stats RI
ir = 2
dist=dists[ir]
srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata(statsdata, dist=dist, srctype=srctype, ptype=ptype, RorI=2)
fig.savefig(f"{srctype}_{ptype}_RI.svg", bbox_inches='tight')
# END plot stats RI
# -------------------------------------------------------------------



# -------------------------------------------------------------------
# BEGIN DEPSRC 0.0 DGRN
depsrc = 0.0
deprcv = 0.0
pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")

pymod.compute_grn(depsrc=depsrc, deprcv=deprcv, dists=dists, nt=500, dt=0.02, converg_method='none', statsidxs=[50,100])

statsdir = f"GRN_grtstats/milrow_{format_float(depsrc)}_{format_float(deprcv)}"
statsdata = pygrt.utils.read_statsfile(f"{statsdir}/K_0050_*")

dist=10
srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata(statsdata, dist=dist, srctype=srctype, ptype=ptype, RorI=2)
fig.savefig(f"{srctype}_{ptype}_{depsrc}_RI.svg", bbox_inches='tight')
# END DEPSRC 0.0 DGRN
# -------------------------------------------------------------------

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in ["GRN", "GRN_grtstats"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
