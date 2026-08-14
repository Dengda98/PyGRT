import shutil
from pathlib import Path

import pygrt

dist = 10.0
depsrc = 2.0
deprcv = 3.0
nt = 600
dt = 0.02
modname = "../milrow"

pymod = pygrt.PyModel1D(grn="GRN", modelpath=modname)

pymod.compute_grn(depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt)
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt, calc_upar=True,
)
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt,
    zeta=0.6, upsampling_n=10,
)
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt,
    freqband=[1, 10],
)
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt, Length=20,
)
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt,
    k0=4, ampk=1.2, keps=1e-3, vmin_ref=1.5,
)

pymod.compute_grn(
    depsrc=depsrc, deprcv=0.0, distarr=2000, nt=1400, dt=1.0, safilonTol=1e-3,
)
pymod.compute_grn(
    depsrc=depsrc, deprcv=0.0, distarr=2000, nt=1400, dt=1.0,
    safilonTol=1e-3, keepAllFreq=True,
)

pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt,
    Length=20, converg_method="DCM",
)
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt,
    Length=20, converg_method="PTAM",
)
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt,
    Length=20, converg_method="none",
)

pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt,
    Length=20, statsidxs=[1, 10, 20],
)

# multi distances
pymod.compute_grn(
    depsrc=depsrc, deprcv=0.0, distarr=[6, 8, 10], nt=nt, dt=dt,
)

# boundary condition
pymod = pygrt.PyModel1D(grn="GRN", modelpath=modname, topbound="free", botbound="free")
pymod.compute_grn(depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt)

pymod = pygrt.PyModel1D(grn="GRN", modelpath=modname, topbound="halfspace", botbound="free")
pymod.compute_grn(depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt)

pymod = pygrt.PyModel1D(grn="GRN", modelpath=modname, topbound="rigid", botbound="rigid")
pymod.compute_grn(depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt)

for name in ["GRN", "GRN_grtstats"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
