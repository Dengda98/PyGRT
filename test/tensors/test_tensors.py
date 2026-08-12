import shutil
from pathlib import Path

import pygrt

dist = 10.0
depsrc = 2.0
deprcv = 3.0
nt = 600
dt = 0.02
modname = "../milrow"
az = 22.0

pymod = pygrt.PyModel1D(modname)
pymod.set_dynamic_grn_path("GRN")
pymod.compute_grn(
    depsrc=depsrc, deprcv=deprcv, distarr=dist, nt=nt, dt=dt, calc_upar=True,
)

pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn", source="EX",
    calc_upar=True,
)
pygrt.utils.compute_strain("syn")
pygrt.utils.compute_stress("syn")
pygrt.utils.compute_rotation("syn")

pymod.compute_syn(
    dist=dist, azimuth=az, scale=1e20, output_path="syn_zne", source="EX",
    zne=True, calc_upar=True,
)
pygrt.utils.compute_strain("syn_zne")
pygrt.utils.compute_stress("syn_zne")
pygrt.utils.compute_rotation("syn_zne")

for name in ["GRN", "syn", "syn_zne"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
