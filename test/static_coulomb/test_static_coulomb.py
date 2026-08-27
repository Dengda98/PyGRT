import numpy as np
import pygrt
from scipy.io import netcdf_file


def read_nc(name):
    with netcdf_file(name, mode="r", mmap=False) as nc:
        data = {key: np.array(value[:], dtype=float, copy=True) for key, value in nc.variables.items()}
        dimensions = set(nc.dimensions)
    return data, dimensions


def check_coulomb(name, friction, dimensions):
    data, actual_dimensions = read_nc(name)
    assert actual_dimensions == dimensions
    assert data["coulomb"].shape == data["sigma_n"].shape
    assert data["coulomb"].shape == data["tau_s"].shape
    assert np.all(np.isfinite(data["coulomb"]))
    np.testing.assert_allclose(
        data["coulomb"], data["tau_s"] + friction * data["sigma_n"],
        rtol=1e-12,
        atol=1e-8,
    )


# 通过 Python API 重新执行一次静态投影和库伦应力写入
pygrt.utils.static_sproj("grid_zne.nc", strike=33.0, dip=44.0, rake=55.0)
pygrt.utils.static_coulomb("grid_zne.nc", 0.6)

check_coulomb("grid_zrt.nc", 0.8, {"north", "east"})
check_coulomb("grid_zne.nc", 0.6, {"north", "east"})
check_coulomb("points_zrt.nc", 0.25, {"point"})
check_coulomb("finite_zrt.nc", 0.4, {"point", "nfault"})

print("test_coulomb.py: all checks passed")
