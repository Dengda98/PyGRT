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
# [('k', '<f8'), ('EX_q', '<c16'), ('EX_w', '<c16'), ('VF_q', '<c16'), ('VF_w', '<c16'), ('HF_q', '<c16'), ('HF_w', '<c16'), ('HF_v', '<c16'), ('DD_q', '<c16'), ('DD_w', '<c16'), ('DS_q', '<c16'), ('DS_w', '<c16'), ('DS_v', '<c16'), ('SS_q', '<c16'), ('SS_w', '<c16'), ('SS_v', '<c16')]
# END read statsfile
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# BEGIN plot stats
ir = 2
dist=distarr[ir]
srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata(statsdata, dist=dist, srctype=srctype, ptype=ptype)
fig.savefig(f"{srctype}_{ptype}.svg", bbox_inches='tight')
# END plot stats
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# BEGIN plot stats RI
ir = 2
dist=distarr[ir]
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
pymod = pygrt.PyModel1D(modarr, depsrc=depsrc, deprcv=deprcv)

stgrnLst = pymod.compute_grn(
    distarr=distarr, nt=500, dt=0.02, converg_method='none',
    statsfile=f"pygrtstats_{depsrc}_{deprcv}", statsidxs=[50,100]
)

statsdata = pygrt.utils.read_statsfile(f"pygrtstats_{depsrc}_{deprcv}/K_0050_*")

dist=10
srctype="SS"
ptype="0"
fig, ax = pygrt.utils.plot_statsdata(statsdata, dist=dist, srctype=srctype, ptype=ptype, RorI=2)
fig.savefig(f"{srctype}_{ptype}_{depsrc}_RI.svg", bbox_inches='tight')
# END DEPSRC 0.0 DGRN
# -------------------------------------------------------------------

