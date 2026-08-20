"""Compare PyGRT Okada finite-fault results with fixed Coulomb references."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat, netcdf_file


DISPLACEMENT_FIELDS = ("Z", "N", "E")
DERIVATIVE_FIELDS = ("zZ", "zN", "zE", "nZ", "nN", "nE", "eZ", "eN", "eE")
ALL_FIELDS = DISPLACEMENT_FIELDS + DERIVATIVE_FIELDS
LABELS = ("kode_100", "kode_200", "kode_300", "kode_400", "kode_500")

# Coulomb input_open reads ELEMENT with %f32, whereas PyGRT reads the decimal
# finite-fault values in double precision.  This tolerance covers that input
# precision difference while still detecting an implementation error.
COULOMB_INPUT_RTOL = 2.0e-4
COULOMB_INPUT_ATOL = 2.0e-12


def read_reference(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    data = loadmat(path)
    norths = np.asarray(data["norths"], dtype=np.float64).reshape(-1)
    easts = np.asarray(data["easts"], dtype=np.float64).reshape(-1)
    depths = np.asarray(data["depths"], dtype=np.float64).reshape(-1)
    if depths.size == 1:
        depths = np.full(norths.size, depths[0], dtype=np.float64)
    fields = {
        name: np.asarray(data[name], dtype=np.float64).squeeze().reshape(-1, order="F")
        for name in ALL_FIELDS
    }
    return norths, easts, depths, fields


def read_pygrt(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    with netcdf_file(path, mmap=False) as nc:
        norths = np.array(nc.variables["north"].data, dtype=np.float64, copy=True).reshape(-1)
        easts = np.array(nc.variables["east"].data, dtype=np.float64, copy=True).reshape(-1)
        depths = np.array(nc.variables["depth"].data, dtype=np.float64, copy=True).reshape(-1)
        fields = {
            name: np.array(nc.variables[name].data, dtype=np.float64, copy=True).reshape(-1)
            for name in ALL_FIELDS
        }
    return norths, easts, depths, fields


def compare_case(label: str, reference_root: Path, result_root: Path) -> None:
    norths, easts, depths, expected = read_reference(
        reference_root / "result" / label / "coulomb.mat"
    )
    actual_norths, actual_easts, actual_depths, actual = read_pygrt(
        result_root / label / "pygrt.nc"
    )

    expected_norths = np.tile(norths, easts.size)
    expected_easts = np.repeat(easts, norths.size)
    np.testing.assert_allclose(actual_norths, expected_norths, rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(actual_easts, expected_easts, rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(actual_depths, depths[0], rtol=0.0, atol=1.0e-12)

    max_abs = 0.0
    max_relative = 0.0
    max_field = ""
    for name in ALL_FIELDS:
        difference = np.abs(actual[name] - expected[name])
        error = float(np.max(difference))
        scale = max(float(np.max(np.abs(expected[name]))), 1.0e-300)
        relative = error / scale
        tolerance = COULOMB_INPUT_RTOL * scale + COULOMB_INPUT_ATOL
        if error > tolerance:
            raise AssertionError(
                f"{label} {name}: max_abs={error:.9g}, max_relative={relative:.9g}, "
                f"allowed={tolerance:.9g}"
            )
        if error > max_abs:
            max_abs = error
            max_field = name
        max_relative = max(max_relative, relative)

    print(
        f"{label}: max_abs={max_abs:.9g} ({max_field}), "
        f"max_relative={max_relative:.9g}"
    )


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: test_okada_coulomb.py <reference-root> <result-root>")

    reference_root = Path(sys.argv[1])
    result_root = Path(sys.argv[2])
    for label in LABELS:
        compare_case(label, reference_root, result_root)
    print("test_okada_coulomb.py: all Coulomb/PyGRT checks passed")


if __name__ == "__main__":
    main()
