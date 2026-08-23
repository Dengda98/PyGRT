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


def read_fault_receiver(path):
    with netcdf_file(path, mmap=False) as dataset:
        layout = getattr(dataset, "layout")
        if isinstance(layout, bytes):
            layout = layout.decode()
        assert layout == "points"
        point = dataset.dimensions["point"]
        nfault = dataset.dimensions["nfault"]
        offset = np.array(dataset.variables["offset"].data, dtype=np.int64, copy=True)
        assert offset.shape == (nfault,)
        assert offset[-1] == point
        point_fields = {
            name: np.array(variable.data, dtype=np.float64, copy=True)
            for name, variable in dataset.variables.items()
            if variable.data.ndim == 1 and variable.data.shape[0] == point
        }
        faults = []
        for ifault in range(nfault):
            start = 0 if ifault == 0 else offset[ifault - 1]
            end = offset[ifault]
            fields = {name: values[start:end] for name, values in point_fields.items()}
            faults.append({
                "coordinates": np.column_stack([
                    fields["north"], fields["east"], fields["depth"],
                ]),
                "fields": fields,
            })
        geometry = {
            name: np.array(dataset.variables[name].data, dtype=np.float64, copy=True)
            for name in ("strike", "dip", "rake", "offset", "stksize", "dipsize")
        }
    return faults, geometry


def write_receiver_points(path, faults):
    coordinates = np.concatenate([fault["coordinates"] for fault in faults])
    lines = ["# north east depth\n"]
    lines.extend(f"{north:.15g} {east:.15g} {depth:.15g}\n" for north, east, depth in coordinates)
    path.write_text("".join(lines))


def assert_faults_match_points(fault_path, points_path, field_names, expected_subsizes=(9, 12)):
    faults, geometry = read_fault_receiver(fault_path)
    assert [fault["coordinates"].shape[0] for fault in faults] == list(expected_subsizes)
    np.testing.assert_array_equal(geometry["offset"], np.cumsum(expected_subsizes))
    assert geometry["rake"][1] == -999.0
    with netcdf_file(points_path, mmap=False) as dataset:
        coordinates = np.column_stack([
            np.array(dataset.variables[name].data, dtype=np.float64, copy=True)
            for name in ("north", "east", "depth")
        ])
        expected_coordinates = np.concatenate([fault["coordinates"] for fault in faults])
        np.testing.assert_allclose(coordinates, expected_coordinates)
        for name in field_names:
            expected = np.concatenate([fault["fields"][name] for fault in faults])
            np.testing.assert_allclose(
                expected, np.array(dataset.variables[name].data, dtype=np.float64, copy=True),
                rtol=1.0e-12, atol=1.0e-12,
            )


def main():
    # The shell commands above provide direct CLI coverage for both output modes
    assert_zne_zrt_equivalent("okada_ff_zne_cli.nc", "okada_ff_zrt_cli.nc")
    assert_receiver_geometry("okada_q6.nc")

    zne_field_names = (
        "Z", "N", "E", "zZ", "zN", "zE", "nZ", "nN", "nE", "eZ", "eN", "eE",
    )
    faults, fault_geometry = read_fault_receiver("okada_rf.nc")
    np.testing.assert_array_equal(fault_geometry["stksize"], [3, 3])
    np.testing.assert_array_equal(fault_geometry["dipsize"], [3, 4])
    rcv_fault_q = Path("rcv_faults_q.txt")
    write_receiver_points(rcv_fault_q, faults)
    compute_okada(
        modelparams=(6.0, 3.464, 2.7),
        depsrc=10.0,
        recv_points=rcv_fault_q,
        output_path="okada_rq.nc",
        scale=1.0e12,
        scale_with_mu=True,
        zne=True,
        calc_upar=True,
    )
    assert_faults_match_points("okada_rf.nc", "okada_rq.nc", zne_field_names)
    default_faults, default_geometry = read_fault_receiver("okada_rf_default.nc")
    assert [fault["coordinates"].shape[0] for fault in default_faults] == [1, 1]
    np.testing.assert_array_equal(default_geometry["offset"], [1, 2])
    np.testing.assert_array_equal(default_geometry["stksize"], [1, 1])
    np.testing.assert_array_equal(default_geometry["dipsize"], [1, 1])
    assert default_geometry["rake"][1] == -999.0

    compute_okada(
        modelparams=(6.0, 3.464, 2.7),
        depsrc=10.0,
        rcv_fault="rcv_faults.inp",
        rcv_fault_size=(0.75, 0.75),
        output_path="okada_rf_api.nc",
        scale=1.0e12,
        scale_with_mu=True,
        zne=True,
        calc_upar=True,
    )
    api_faults, _ = read_fault_receiver("okada_rf_api.nc")
    for shell_fault, api_fault in zip(faults, api_faults):
        np.testing.assert_allclose(shell_fault["coordinates"], api_fault["coordinates"])
        for name in zne_field_names:
            np.testing.assert_allclose(shell_fault["fields"][name], api_fault["fields"][name])

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
        src_fault="cfaults.inp",
        calc_upar=True,
    )
    compute_okada(
        modelparams=(6.0, 3.464, 2.7),
        deprcv=0.0,
        norths=(-5.0, 5.0, 0.5),
        easts=(-5.0, 5.0, 0.5),
        output_path="okada_python_ff_zne.nc",
        src_fault="cfaults.inp",
        zne=True,
        calc_upar=True,
    )
    assert_zne_zrt_equivalent("okada_python_ff_zne.nc", "okada_python_ff_zrt.nc")

    print("test_okada.py: all checks passed")

    Path("okada_python.nc").unlink(missing_ok=True)
    Path("okada_python_ff.nc").unlink(missing_ok=True)
    Path("okada_python_ff_zrt.nc").unlink(missing_ok=True)
    Path("okada_python_ff_zne.nc").unlink(missing_ok=True)
    Path("okada_rq.nc").unlink(missing_ok=True)
    Path("okada_rf_default.nc").unlink(missing_ok=True)
    Path("okada_rf_api.nc").unlink(missing_ok=True)
    Path("rcv_faults_q.txt").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
