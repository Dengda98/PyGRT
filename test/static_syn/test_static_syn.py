from pathlib import Path

import pygrt

depsrc = 2.0
deprcv = 0.0
norths = [-3.0, 3.0, 0.2]
easts = [-2.0, 2.0, 0.2]
modname = "../milrow"

pymod = pygrt.PyModel1D(modname)
pymod.set_static_grn_path("stgrn.nc")
pymod.compute_static_grn(
    depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, calc_upar=True,
)

pymod.compute_static_syn(
    scale=1e20, output_path="stsyn.nc", source="EX",
)
pymod.compute_static_syn(
    scale=1e16, output_path="stsyn.nc", source="SF", force=(-1, 2, -4),
)
pymod.compute_static_syn(
    scale=1e20, output_path="stsyn.nc", source="DC",
    strike=33, dip=44, rake=55,
)
pymod.compute_static_syn(
    scale=1e20, output_path="stsyn.nc", source="TS", strike=33, dip=44,
)
pymod.compute_static_syn(
    scale=1e20, output_path="stsyn.nc", source="MT",
    moment_tensor=(1, -2, -5, 0.5, 3, 1.2),
)

for name in ["stgrn.nc", "stsyn.nc"]:
    p = Path(name)
    if p.is_file():
        p.unlink(missing_ok=True)
