import shutil
from pathlib import Path

import pygrt

depsrc = 2.0
deprcv = 0.0
# norths/easts 各为 start/stop/step (km)
norths = [-3.0, 3.0, 0.2]
easts = [-2.0, 2.0, 0.2]
modname = "../milrow"

pymod = pygrt.PyModel1D(modname)
pymod.set_static_grn_path("stgrn.nc")

pymod.compute_static_grn(depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, calc_upar=True,
)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, Length=20,
)

pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts,
    Length=20, converg_method="DCM",
)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts,
    Length=20, converg_method="PTAM",
)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts,
    Length=20, converg_method="none",
)

pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, k0=4, keps=1e-3,
)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, k0=4, stats=True,
)

# boundary condition
pymod = pygrt.PyModel1D(modname, topbound="free", botbound="free")
pymod.set_static_grn_path("stgrn.nc")
pymod.compute_static_grn(depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)

pymod = pygrt.PyModel1D(modname, topbound="halfspace", botbound="free")
pymod.set_static_grn_path("stgrn.nc")
pymod.compute_static_grn(depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)

pymod = pygrt.PyModel1D(modname, topbound="rigid", botbound="rigid")
pymod.set_static_grn_path("stgrn.nc")
pymod.compute_static_grn(depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts)

for name in ["stgrn.nc", "stgrtstats"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
