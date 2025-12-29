# -------------------------------------------------------------------
# BEGIN DEPSRC 0.0 DGRN
import numpy as np
import pygrt 

modarr = np.loadtxt("milrow")

depsrc = 0.0
deprcv = 0.0
pymod = pygrt.PyModel1D(modarr, depsrc=depsrc, deprcv=deprcv)

distarr = [5,8,10]
# 设置 converg_method='PTAM' 进行收敛
stgrnLst = pymod.compute_grn(
    distarr=distarr, nt=500, dt=0.02, converg_method='PTAM',
    statsfile=f"pygrtstats_{depsrc}_{deprcv}", statsidxs=[50,100]
)
# END DEPSRC 0.0 DGRN
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# BEGIN plot ptam
ir = 2
statsdata1, statsdata2, ptamdata, dist = pygrt.utils.read_statsfile_ptam(f"pygrtstats_{depsrc}_{deprcv}/PTAM_{ir:04d}_*/PTAM_0050_*")

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

modarr = np.loadtxt("milrow")

depsrc = 0.1
deprcv = 0.0
pymod = pygrt.PyModel1D(modarr, depsrc=depsrc, deprcv=deprcv)

xarr = np.array([2.0])
yarr = np.array([2.0])
static_grn = pymod.compute_static_grn(xarr, yarr, converg_method='PTAM', statsfile=f"static_pygrtstats_{depsrc}_{deprcv}")

ir = 0
statsdata1, statsdata2, ptamdata, dist = pygrt.utils.read_statsfile_ptam(f"static_pygrtstats_{depsrc}_{deprcv}/PTAM_{ir:04d}_*/PTAM")

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
