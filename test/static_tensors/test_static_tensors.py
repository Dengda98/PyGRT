from pathlib import Path

import pygrt

depsrc = 2.0
deprcv = 0.0
norths = [-3.0, 3.0, 0.2]
easts = [-2.0, 2.0, 0.2]
modname = "../milrow"

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=modname)
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, calc_upar=True,
)

pymod.compute_static_syn(
    scale=1e20, output_path="stsyn.nc", source="EX", calc_upar=True,
)
pygrt.utils.compute_strain("stsyn.nc")
pygrt.utils.compute_stress("stsyn.nc")
pygrt.utils.compute_rotation("stsyn.nc")

pymod.compute_static_syn(
    scale=1e20, output_path="stsyn_zne.nc", source="EX",
    zne=True, calc_upar=True,
)
pygrt.utils.compute_strain("stsyn_zne.nc")
pygrt.utils.compute_stress("stsyn_zne.nc")
pygrt.utils.compute_rotation("stsyn_zne.nc")

for name in ["stgrn.nc", "stsyn.nc", "stsyn_zne.nc"]:
    p = Path(name)
    if p.is_file():
        p.unlink(missing_ok=True)
