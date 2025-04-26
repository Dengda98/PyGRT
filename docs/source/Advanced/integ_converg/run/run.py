# -------------------------------------------------------------------
# BEGIN DGRN
import numpy as np
import pygrt 

modarr = np.loadtxt("milrow")

depsrc = 2.0
deprcv = 0.0

pymod = pygrt.PyModel1D(modarr, depsrc=depsrc, deprcv=deprcv)

# 通过statsfile参数自定义核函数文件的输出目录, statsidxs指定想输出的频率索引值
distarr = [5,8,10]
stgrnLst = pymod.compute_grn(
    distarr=distarr, nt=500, dt=0.02,
    statsfile=f"pygrtstats_{depsrc}_{deprcv}", statsidxs=[50,100]
)
# END DGRN
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# BEGIN read statsfile
# 可使用通配符简化输入，因为对应索引值下只会有一个文件
# 返回的是自定义类型的numpy数组
statsdata = pygrt.utils.read_statsfile(f"pygrtstats_{depsrc}_{deprcv}/K_0050_*")
print(statsdata.dtype)
# [('k', '<f8'), ('EXP_q0', '<c16'), ('EXP_w0', '<c16'), ('VF_q0', '<c16'), ('VF_w0', '<c16'), ('HF_q1', '<c16'), ('HF_w1', '<c16'), ('HF_v1', '<c16'), ('DC_q0', '<c16'), ('DC_w0', '<c16'), ('DC_q1', '<c16'), ('DC_w1', '<c16'), ('DC_v1', '<c16'), ('DC_q2', '<c16'), ('DC_w2', '<c16'), ('DC_v2', '<c16')]
# END read statsfile
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# BEGIN plot stats
ir = 2
dist=distarr[ir]
srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata(statsdata, dist=dist, srctype=srctype, ptype=ptype)
fig.tight_layout()
fig.savefig(f"{srctype}_{ptype}.png", dpi=100)
# END plot stats
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# BEGIN plot stats RI
ir = 2
dist=distarr[ir]
srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata(statsdata, dist=dist, srctype=srctype, ptype=ptype, RorI=2)
fig.tight_layout()
fig.savefig(f"{srctype}_{ptype}_RI.png", dpi=100)
# END plot stats RI
# -------------------------------------------------------------------



# -------------------------------------------------------------------
# BEGIN DEPSRC 0.0 DGRN
depsrc = 0.0
deprcv = 0.0
pymod = pygrt.PyModel1D(modarr, depsrc=depsrc, deprcv=deprcv)

stgrnLst = pymod.compute_grn(
    distarr=distarr, nt=500, dt=0.02,
    statsfile=f"pygrtstats_{depsrc}_{deprcv}", statsidxs=[50,100]
)

statsdata = pygrt.utils.read_statsfile(f"pygrtstats_{depsrc}_{deprcv}/K_0050_*")

dist=10
srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata(statsdata, dist=dist, srctype=srctype, ptype=ptype, RorI=2)
fig.tight_layout()
fig.savefig(f"{srctype}_{ptype}_{depsrc}_RI.png", dpi=100)
# END DEPSRC 0.0 DGRN
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# BEGIN plot ptam
ir = 2
statsdata1, statsdata2, ptamdata, dist = pygrt.utils.read_statsfile_ptam(f"pygrtstats_{depsrc}_{deprcv}/PTAM_{ir:04d}_*/PTAM_0050_*")

srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata_ptam(statsdata1, statsdata2, ptamdata, dist=dist, srctype=srctype, ptype=ptype, RorI=2)
fig.tight_layout()
fig.savefig(f"{srctype}_{ptype}_{depsrc}_ptam_RI.png", dpi=100)
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
static_grn = pymod.compute_static_grn(xarr, yarr, statsfile=f"static_pygrtstats_{depsrc}_{deprcv}")

ir = 0
statsdata1, statsdata2, ptamdata, dist = pygrt.utils.read_statsfile_ptam(f"static_pygrtstats_{depsrc}_{deprcv}/PTAM_{ir:04d}_*/PTAM")

srctype="SS"
ptype="0"
# 只使用离散波数积分的积分变化
fig, ax = pygrt.utils.plot_statsdata(statsdata1, dist=dist, srctype=srctype, ptype=ptype, RorI=True)
fig.tight_layout()
fig.savefig(f"{srctype}_{ptype}_{depsrc}_static.png", dpi=100)

# 使用了峰谷平均法的积分变化
fig, ax = pygrt.utils.plot_statsdata_ptam(statsdata1, statsdata2, ptamdata, dist=dist, srctype=srctype, ptype=ptype, RorI=True)
fig.tight_layout()
fig.savefig(f"{srctype}_{ptype}_{depsrc}_ptam_static.png", dpi=100)

# END SGRN
# -------------------------------------------------------------------
