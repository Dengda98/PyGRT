# -------------------------------------------------------------------
# BEGIN DEPSRC 0.0 DGRN
import numpy as np
import pygrt 
from pygrt.cli import format_float

depsrc = 0.0
deprcv = 0.0
pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")

distarr = [5,8,10]
# 设置 converg_method='PTAM' 进行收敛
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv,
    distarr=distarr, nt=500, dt=0.02, converg_method='PTAM', k0=2, ampk=1.2, use_kmax_ref=True,
    statsidxs=[50,100],
)
# END DEPSRC 0.0 DGRN
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# BEGIN plot ptam
ir = 2
statsdata1, statsdata2, ptamdata, dist = pygrt.utils.read_statsfile_ptam(
    f"GRN_grtstats/milrow_{format_float(depsrc)}_{format_float(deprcv)}/PTAM_{ir:04d}_*/PTAM_0050_*"
)

srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata_ptam(statsdata1, statsdata2, ptamdata, dist=dist, srctype=srctype, ptype=ptype, RorI=2)
fig.savefig(f"{srctype}_{ptype}_{depsrc}_ptam_RI.svg", bbox_inches='tight')
# END plot ptam
# -------------------------------------------------------------------



# -------------------------------------------------------------------
# BEGIN SGRN
import numpy as np
import pygrt
from pygrt.cli import format_float

depsrc = 0.05
deprcv = 0.0
pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath="milrow")

norths = [2.0, 2.0, 1.0]
easts = [2.0, 2.0, 1.0]
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv,
    norths=norths, easts=easts,
    converg_method='PTAM', stats=True, k0=3, use_kmax_ref=True,
)

ir = 0
statsdata1, statsdata2, ptamdata, dist = pygrt.utils.read_statsfile_ptam(
    f"stgrtstats/milrow_{format_float(depsrc)}_{format_float(deprcv)}/PTAM_{ir:04d}_*/PTAM"
)

srctype="SS"
ptype="0"
# 只使用离散波数积分的积分变化
fig, ax = pygrt.utils.plot_statsdata(statsdata1, dist=dist, srctype=srctype, ptype=ptype, RorI=True)
fig.tight_layout()
fig.savefig(f"{srctype}_{ptype}_{depsrc}_static.svg", bbox_inches='tight')

# 使用了峰谷平均法的积分变化
fig, ax = pygrt.utils.plot_statsdata_ptam(statsdata1, statsdata2, ptamdata, dist=dist, srctype=srctype, ptype=ptype, RorI=True)
fig.tight_layout()
fig.savefig(f"{srctype}_{ptype}_{depsrc}_ptam_static.svg", bbox_inches='tight')

# END SGRN
# -------------------------------------------------------------------

# 删除中间计算结果，仅保留成图
import shutil
from pathlib import Path
for name in ["GRN", "GRN_grtstats", "stgrn.nc", "stgrtstats"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
