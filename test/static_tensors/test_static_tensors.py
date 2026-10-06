from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import pygrt
from scipy.io import netcdf_file

depsrc = 2.0
deprcv = 0.0
norths = [-3.0, 3.0, 0.2]
easts = [-2.0, 2.0, 0.2]
modname = "../milrow"

pymod = pygrt.PyModel1D(stgrn="stgrn.nc", modelpath=modname)
pymod.static_greenfn(depsrc=depsrc, deprcv=deprcv, norths=norths, easts=easts, calc_upar=True)

pymod.static_syn(scale=1e20, output_path="stsyn.nc", calc_upar=True)
pygrt.utils.static_strain("stsyn.nc")
pygrt.utils.static_stress("stsyn.nc")
pygrt.utils.static_rotation("stsyn.nc")

pymod.static_syn(scale=1e20, output_path="stsyn_zne.nc", zne=True, calc_upar=True)
pygrt.utils.static_strain("stsyn_zne.nc")
pygrt.utils.static_stress("stsyn_zne.nc")
pygrt.utils.static_rotation("stsyn_zne.nc")

# 坐标维度与布局不一致时，在读取数组和写入后处理结果之前拒绝文件
with TemporaryDirectory(prefix="pygrt-static-nc-") as temporary:
    for layout, coordinate_dims, expected_error in [
        ("points", ("point", "extra"), "must have one dimension"),
        ("points", ("extra",), "has an unexpected dimension"),
        ("points", (), "must have one dimension"),
        ("grid", ("east",), "has an unexpected dimension"),
    ]:
        path = Path(temporary)/"invalid_coordinates.nc"
        with netcdf_file(path, "w") as dataset:
            dataset.layout = layout
            dataset.rot2ZNE = 1
            dataset.calc_upar = 1
            for name, count in [("point", 2), ("extra", 3), ("north", 2), ("east", 3)]:
                dataset.createDimension(name, count)
            north = dataset.createVariable("north", "d", coordinate_dims)
            north.data[...] = 0
            east = dataset.createVariable("east", "d", ("east",) if layout == "grid" else ("point",))
            east[:] = 0
        before = path.read_bytes()
        result = subprocess.run(["grt", "static_strain", str(path)], capture_output=True, text=True)
        assert result.returncode != 0
        assert expected_error in result.stdout+result.stderr
        assert path.read_bytes() == before

    path = Path(temporary)/"empty_points.nc"
    with netcdf_file(path, "w") as dataset:
        dataset.layout = "points"
        dataset.rot2ZNE = 1
        dataset.calc_upar = 1
        dataset.createDimension("point", 0)
    before = path.read_bytes()
    result = subprocess.run(["grt", "static_strain", str(path)], capture_output=True, text=True)
    assert result.returncode != 0
    assert "dimension point is empty" in result.stdout+result.stderr
    assert path.read_bytes() == before

for name in ["stgrn.nc", "stsyn.nc", "stsyn_zne.nc"]:
    p = Path(name)
    if p.is_file():
        p.unlink(missing_ok=True)
