from pathlib import Path

import numpy as np
from pygrt.utils import compute_okada
from scipy.io import netcdf_file


def read_okada_fields(path):
    with netcdf_file(path, mmap=False) as dataset:
        variables = {
            name: np.array(variable.data, dtype=np.float64, copy=True)
            for name, variable in dataset.variables.items()
        }
        attributes = {
            name: getattr(dataset, name)
            for name in ("rot2ZNE", "calc_upar")
            if hasattr(dataset, name)
        }
    return variables, attributes


def assert_zne_zrt_equivalent(zne_path, zrt_path):
    zne, zne_attributes = read_okada_fields(zne_path)
    zrt, zrt_attributes = read_okada_fields(zrt_path)
    assert zne_attributes["rot2ZNE"] == 1
    assert zrt_attributes["rot2ZNE"] == 0
    assert zne_attributes["calc_upar"] == 1
    assert zrt_attributes["calc_upar"] == 1
    np.testing.assert_allclose(zne["north"], zrt["north"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(zne["east"], zrt["east"], rtol=0.0, atol=0.0)

    north = zrt["north"][:, None]
    east = zrt["east"][None, :]
    distance = np.hypot(north, east)
    theta = np.where(distance == 0.0, 0.0, np.arctan2(east, north))
    cosine = np.cos(theta)
    sine = np.sin(theta)

    np.testing.assert_allclose(zne["Z"], zrt["Z"], rtol=1.0e-12, atol=1.0e-12)
    np.testing.assert_allclose(
        zne["N"], zrt["R"] * cosine - zrt["T"] * sine,
        rtol=1.0e-12, atol=1.0e-12,
    )
    np.testing.assert_allclose(
        zne["E"], zrt["R"] * sine + zrt["T"] * cosine,
        rtol=1.0e-12, atol=1.0e-12,
    )

    r = distance * 1.0e5
    s00, s01, s02 = zrt["zZ"], zrt["zR"], zrt["zT"]
    s10, s11, s12 = zrt["rZ"], zrt["rR"], zrt["rT"]
    s20, s21, s22 = zrt["tZ"], zrt["tR"], zrt["tT"]
    ur_over_r = np.array(s11, copy=True)
    ut_over_r = np.array(s12, copy=True)
    nonzero = r != 0.0
    ur_over_r[nonzero] = zrt["R"][nonzero] / r[nonzero]
    ut_over_r[nonzero] = zrt["T"][nonzero] / r[nonzero]
    converted = {
        "zZ": s00,
        "zN": s01 * cosine - s02 * sine,
        "zE": s01 * sine + s02 * cosine,
        "nZ": s10 * cosine - s20 * sine,
        "eZ": s10 * sine + s20 * cosine,
        "nN": s11 * cosine**2 + s22 * sine**2 - (s12 + s21) * sine * cosine
        + ur_over_r * sine**2 + ut_over_r * sine * cosine,
        "nE": s12 * cosine**2 - s21 * sine**2 + (s11 - s22) * sine * cosine
        - ur_over_r * sine * cosine + ut_over_r * sine**2,
        "eN": s21 * cosine**2 - s12 * sine**2 + (s11 - s22) * sine * cosine
        - ur_over_r * sine * cosine - ut_over_r * cosine**2,
        "eE": s22 * cosine**2 + s11 * sine**2 + (s12 + s21) * sine * cosine
        + ur_over_r * cosine**2 - ut_over_r * sine * cosine,
    }
    for name, expected in converted.items():
        np.testing.assert_allclose(zne[name], expected, rtol=2.0e-11, atol=1.0e-12)


def assert_receiver_geometry(path):
    fields, _ = read_okada_fields(path)
    for name in ("strike", "dip", "rake"):
        assert name in fields
        assert fields[name].shape == (3,)
    np.testing.assert_allclose(fields["strike"], [10.0, 40.0, 70.0])
    np.testing.assert_allclose(fields["dip"], [20.0, 50.0, 80.0])
    np.testing.assert_allclose(fields["rake"], [30.0, 60.0, 90.0])


def main():
    # The shell commands above provide direct CLI coverage for both output modes
    assert_zne_zrt_equivalent("okada_ff_zne_cli.nc", "okada_ff_zrt_cli.nc")
    assert_receiver_geometry("okada_q6.nc")

    compute_okada(
        modelparams=(6.0, 3.464, 2.7),
        depsrc=50.0,
        deprcv=0.0,
        norths=(-5.0, 5.0, 0.5),
        easts=(-5.0, 5.0, 0.5),
        output_path="okada_python.nc",
        scale=1.0e12,
    )
    compute_okada(
        modelparams=(6.0, 3.464, 2.7),
        deprcv=0.0,
        norths=(-5.0, 5.0, 0.5),
        easts=(-5.0, 5.0, 0.5),
        output_path="okada_python_ff_zrt.nc",
        finite_fault="cfaults.inp",
        calc_upar=True,
    )
    compute_okada(
        modelparams=(6.0, 3.464, 2.7),
        deprcv=0.0,
        norths=(-5.0, 5.0, 0.5),
        easts=(-5.0, 5.0, 0.5),
        output_path="okada_python_ff_zne.nc",
        finite_fault="cfaults.inp",
        zne=True,
        calc_upar=True,
    )
    assert_zne_zrt_equivalent("okada_python_ff_zne.nc", "okada_python_ff_zrt.nc")

    print("test_okada.py: all checks passed")

    Path("okada_python.nc").unlink(missing_ok=True)
    Path("okada_python_ff.nc").unlink(missing_ok=True)
    Path("okada_python_ff_zrt.nc").unlink(missing_ok=True)
    Path("okada_python_ff_zne.nc").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
