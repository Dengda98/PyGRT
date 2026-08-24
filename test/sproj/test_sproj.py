from pathlib import Path

import numpy as np
from scipy.io import netcdf_file


HERE = Path(__file__).resolve().parent


def read_nc(name):
    with netcdf_file(HERE / name, mode="r", mmap=False) as nc:
        data = {key: np.array(value[:], dtype=float, copy=True) for key, value in nc.variables.items()}
        dimensions = set(nc.dimensions)
    return data, dimensions


def receiver_vectors(strike, dip, rake):
    strike = np.deg2rad(strike)
    dip = np.deg2rad(dip)
    rake = np.deg2rad(rake)
    nvec = np.array(
        [-np.sin(strike) * np.sin(dip), np.cos(strike) * np.sin(dip), np.cos(dip)]
    )
    tvec = np.array(
        [
            np.cos(rake) * np.cos(strike) + np.sin(rake) * np.cos(dip) * np.sin(strike),
            np.cos(rake) * np.sin(strike) - np.sin(rake) * np.cos(dip) * np.cos(strike),
            np.sin(rake) * np.sin(dip),
        ]
    )
    return nvec, tvec


def project_zne(data, strike, dip, rake, index):
    stress = np.array(
        [
            data["stress_ZZ"].flat[index],
            data["stress_ZN"].flat[index],
            data["stress_ZE"].flat[index],
            data["stress_NN"].flat[index],
            data["stress_NE"].flat[index],
            data["stress_EE"].flat[index],
        ]
    )
    nvec, tvec = receiver_vectors(strike, dip, rake)
    traction = np.array(
        [
            stress[3] * nvec[0] + stress[4] * nvec[1] + stress[1] * nvec[2],
            stress[4] * nvec[0] + stress[5] * nvec[1] + stress[2] * nvec[2],
            stress[1] * nvec[0] + stress[2] * nvec[1] + stress[0] * nvec[2],
        ]
    )
    return np.dot(traction, nvec), np.dot(traction, tvec)


def project_zrt(data, strike, dip, rake, index):
    north = np.repeat(data["north"], data["east"].size).flat[index]
    east = np.tile(data["east"], data["north"].size).flat[index]
    theta = np.arctan2(east, north) if np.hypot(north, east) > 1e-8 else 0.0
    stress = np.array(
        [
            data["stress_ZZ"].flat[index],
            data["stress_ZR"].flat[index],
            data["stress_ZT"].flat[index],
            data["stress_RR"].flat[index],
            data["stress_RT"].flat[index],
            data["stress_TT"].flat[index],
        ]
    )
    st = np.sin(-theta)
    ct = np.cos(-theta)
    sst = st * st
    cct = ct * ct
    sct = st * ct
    zne = np.array(
        [
            stress[0],
            stress[1] * ct + stress[2] * st,
            -stress[1] * st + stress[2] * ct,
            stress[3] * cct + stress[5] * sst + 2.0 * stress[4] * sct,
            (stress[5] - stress[3]) * sct + stress[4] * (cct - sst),
            stress[3] * sst + stress[5] * cct - 2.0 * stress[4] * sct,
        ]
    )
    zne_data = {
        "stress_ZZ": np.asarray([zne[0]]),
        "stress_ZN": np.asarray([zne[1]]),
        "stress_ZE": np.asarray([zne[2]]),
        "stress_NN": np.asarray([zne[3]]),
        "stress_NE": np.asarray([zne[4]]),
        "stress_EE": np.asarray([zne[5]]),
    }
    return project_zne(zne_data, strike, dip, rake, 0)


grid_zrt, grid_zrt_dims = read_nc("grid_zrt.nc")
grid_zne, grid_zne_dims = read_nc("grid_zne.nc")
assert grid_zrt_dims == {"north", "east"}
assert grid_zne_dims == {"north", "east"}
assert grid_zrt["sigma_n"].shape == (5, 3)
assert grid_zrt["tau_s"].shape == (5, 3)

expected_sigma, expected_tau = project_zne(grid_zne, 33.0, 44.0, 55.0, 0)
np.testing.assert_allclose(grid_zne["sigma_n"].flat[0], expected_sigma, rtol=1e-12, atol=1e-8)
np.testing.assert_allclose(grid_zne["tau_s"].flat[0], expected_tau, rtol=1e-12, atol=1e-8)
expected_sigma, expected_tau = project_zrt(grid_zrt, 33.0, 44.0, 55.0, 0)
np.testing.assert_allclose(grid_zrt["sigma_n"].flat[0], expected_sigma, rtol=1e-12, atol=1e-8)
np.testing.assert_allclose(grid_zrt["tau_s"].flat[0], expected_tau, rtol=1e-12, atol=1e-8)

plain_manual, plain_dims = read_nc("points_plain_manual.nc")
assert plain_dims == {"point"}
assert plain_manual["sigma_n"].shape == (3,)
assert plain_manual["tau_s"].shape == (3,)

points_from_file, _ = read_nc("points_geometry_from_file.nc")
points_manual, _ = read_nc("points_geometry.nc")
assert points_from_file["strike"].shape == (3,)
assert points_from_file["dip"].shape == (3,)
assert points_from_file["rake"].shape == (3,)
assert not np.allclose(points_from_file["sigma_n"], points_manual["sigma_n"])
assert not np.allclose(points_from_file["tau_s"], points_manual["tau_s"])

points_q, _ = read_nc("points_plain.nc")
assert not np.allclose(plain_manual["sigma_n"], points_q["sigma_n"])
assert not np.allclose(plain_manual["tau_s"], points_q["tau_s"])

finite_undefined, finite_undefined_dims = read_nc("finite_undefined.nc")
finite_undefined_partial, _ = read_nc("finite_undefined_partial.nc")
finite_defined, finite_defined_dims = read_nc("finite_defined.nc")
finite_no_rake, finite_no_rake_dims = read_nc("finite_no_rake.nc")
assert "nfault" in finite_undefined_dims
assert "nfault" in finite_defined_dims
assert "nfault" in finite_no_rake_dims
assert finite_undefined["sigma_n"].shape == finite_undefined["north"].shape
assert finite_undefined_partial["sigma_n"].shape == finite_undefined_partial["north"].shape
assert finite_defined["sigma_n"].shape == finite_defined["north"].shape
assert finite_no_rake["sigma_n"].shape == finite_no_rake["north"].shape
assert np.all(np.isfinite(finite_undefined["sigma_n"]))
assert np.all(np.isfinite(finite_undefined["tau_s"]))
assert not np.allclose(finite_undefined_partial["tau_s"][4:], finite_undefined["tau_s"][4:])
assert np.all(np.isfinite(finite_defined["sigma_n"]))
assert np.all(np.isfinite(finite_defined["tau_s"]))
assert np.all(np.isfinite(finite_no_rake["sigma_n"]))
assert np.all(np.isfinite(finite_no_rake["tau_s"]))

print("test_sproj.py: all checks passed")
