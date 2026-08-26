from pathlib import Path
import re

import numpy as np
from scipy.io import netcdf_file

import pygrt.utils


EARTH_RADIUS_KM = 6371.0
HERE = Path(__file__).resolve().parent


def local_to_geo(north, east, lon0, lat0):
    km_per_lat_deg = EARTH_RADIUS_KM * np.pi / 180.0
    km_per_lon_deg = EARTH_RADIUS_KM * np.cos(np.deg2rad(lat0)) * np.pi / 180.0
    latitude = lat0 + np.asarray(north) / km_per_lat_deg
    longitude = lon0 + np.asarray(east) / km_per_lon_deg
    longitude = (longitude + 180.0) % 360.0 - 180.0
    return latitude, longitude


def geo_to_local(latitude, longitude, lon0, lat0):
    km_per_lat_deg = EARTH_RADIUS_KM * np.pi / 180.0
    km_per_lon_deg = EARTH_RADIUS_KM * np.cos(np.deg2rad(lat0)) * np.pi / 180.0
    longitude_delta = (np.asarray(longitude) - lon0 + 180.0) % 360.0 - 180.0
    north = (np.asarray(latitude) - lat0) * km_per_lat_deg
    east = longitude_delta * km_per_lon_deg
    return north, east


def scalar_attribute(dataset, name):
    value = getattr(dataset, name)
    if isinstance(value, np.ndarray):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode()
    return value


def assert_reference(dataset, lon0, lat0):
    np.testing.assert_allclose(float(scalar_attribute(dataset, "lon0")), lon0)
    np.testing.assert_allclose(float(scalar_attribute(dataset, "lat0")), lat0)


def assert_nc_roundtrip(name, layout, north, east, lon0, lat0, input_dims):
    with netcdf_file(HERE / f"{name}_local.nc", mode="r", mmap=False) as before:
        assert set(before.dimensions) == set(input_dims)
        np.testing.assert_allclose(before.variables["north"][:], north)
        np.testing.assert_allclose(before.variables["east"][:], east)
        input_field = np.array(before.variables["field"][:], copy=True)
        input_depth = (
            np.array(before.variables["depth"][:], copy=True)
            if "depth" in before.variables
            else None
        )

    expected_lat, expected_lon = local_to_geo(north, east, lon0, lat0)
    expected_geo_dims = {"lat", "lon"} if layout == "grid" else set(input_dims)
    with netcdf_file(HERE / f"{name}_geo.nc", mode="r", mmap=False) as after:
        assert set(after.dimensions) == expected_geo_dims
        assert set(after.variables).isdisjoint({"north", "east"})
        assert {"lat", "lon"}.issubset(after.variables)
        if layout == "grid":
            assert after.variables["lat"].dimensions == ("lat",)
            assert after.variables["lon"].dimensions == ("lon",)
        else:
            assert after.variables["lat"].dimensions == ("point",)
            assert after.variables["lon"].dimensions == ("point",)
        np.testing.assert_allclose(after.variables["lat"][:], expected_lat)
        np.testing.assert_allclose(after.variables["lon"][:], expected_lon)
        np.testing.assert_array_equal(after.variables["field"][:], input_field)
        if input_depth is not None:
            np.testing.assert_array_equal(after.variables["depth"][:], input_depth)
        assert_reference(after, lon0, lat0)
        assert scalar_attribute(after, "layout") == layout
        assert scalar_attribute(after, "computeType") == "test"

    expected_north, expected_east = geo_to_local(expected_lat, expected_lon, lon0, lat0)
    with netcdf_file(HERE / f"{name}_roundtrip.nc", mode="r", mmap=False) as back:
        assert set(back.dimensions) == set(input_dims)
        expected_north_dims = ("north",) if layout == "grid" else ("point",)
        expected_east_dims = ("east",) if layout == "grid" else ("point",)
        assert back.variables["north"].dimensions == expected_north_dims
        assert back.variables["east"].dimensions == expected_east_dims
        np.testing.assert_allclose(back.variables["north"][:], expected_north, atol=1e-12)
        np.testing.assert_allclose(back.variables["east"][:], expected_east, atol=1e-12)
        np.testing.assert_array_equal(back.variables["field"][:], input_field)
        if input_depth is not None:
            np.testing.assert_array_equal(back.variables["depth"][:], input_depth)
        assert_reference(back, lon0, lat0)


assert_nc_roundtrip(
    "grid",
    "grid",
    np.array([-10.0, 0.0, 10.0]),
    np.array([-200.0, 0.0, 200.0]),
    179.5,
    35.0,
    {"north", "east"},
)

assert_nc_roundtrip(
    "points",
    "points",
    np.array([0.0, 10.0, -20.0]),
    np.array([-200.0, 0.0, 200.0]),
    -179.5,
    -20.0,
    {"point"},
)

assert_nc_roundtrip(
    "finite",
    "points",
    np.array([0.0, 10.0, 20.0, 30.0]),
    np.array([-2.0, -1.0, 1.0, 2.0]),
    10.0,
    40.0,
    {"point", "nfault"},
)

with netcdf_file(HERE / "finite_geo.nc", mode="r", mmap=False) as nc:
    assert nc.variables["strike"].dimensions == ("nfault",)
    assert nc.variables["offset"].dimensions == ("nfault",)
    np.testing.assert_array_equal(nc.variables["offset"][:], [2, 4])


def parse_data_line(line):
    match = re.match(
        r"(?P<prefix>\s*)(?P<first>[^\s]+)(?P<separator>\s+)"
        r"(?P<second>[^\s]+)(?P<tail>.*)",
        line,
    )
    assert match is not None, line
    return match


def assert_text_conversion():
    original_lines = (HERE / "local_points.txt").read_text().splitlines(keepends=True)
    geo_lines = (HERE / "geo_points.txt").read_text().splitlines(keepends=True)
    roundtrip_lines = (HERE / "local_points_roundtrip.txt").read_text().splitlines(keepends=True)
    assert geo_lines[0] == "# 35 179.5\n"
    assert roundtrip_lines[0] == "# 35 179.5\n"
    assert len(geo_lines) == len(original_lines) + 1
    assert len(roundtrip_lines) == len(geo_lines) + 1

    north = []
    east = []
    data_indices = []
    for index, line in enumerate(original_lines):
        if not line.strip() or line.lstrip().startswith("#"):
            assert geo_lines[index + 1] == line
            continue
        match = parse_data_line(line)
        north.append(float(match.group("first")))
        east.append(float(match.group("second")))
        data_indices.append(index)
        converted = parse_data_line(geo_lines[index + 1])
        expected_lat, expected_lon = local_to_geo(
            float(match.group("first")), float(match.group("second")), 179.5, 35.0
        )
        np.testing.assert_allclose(float(converted.group("first")), expected_lat)
        np.testing.assert_allclose(float(converted.group("second")), expected_lon)
        assert converted.group("prefix") == match.group("prefix")
        assert converted.group("tail") == match.group("tail")

    expected_north, expected_east = geo_to_local(
        *local_to_geo(np.array(north), np.array(east), 179.5, 35.0), 179.5, 35.0
    )
    for data_number, index in enumerate(data_indices):
        source = parse_data_line(geo_lines[index + 1])
        converted = parse_data_line(roundtrip_lines[index + 2])
        np.testing.assert_allclose(float(converted.group("first")), expected_north[data_number])
        np.testing.assert_allclose(float(converted.group("second")), expected_east[data_number])
        assert converted.group("prefix") == source.group("prefix")
        assert converted.group("tail") == source.group("tail")
    for index, line in enumerate(geo_lines):
        if not line.strip() or line.lstrip().startswith("#"):
            assert roundtrip_lines[index + 1] == line

    # The Python wrappers map -G and -Q to the corresponding commands
    class CapturedRunner:
        def __init__(self):
            self.commands = []

        def __call__(self, command, **kwargs):
            self.commands.append([str(value) for value in command])

    runner = CapturedRunner()
    original_runner = pygrt.utils.run_grt
    pygrt.utils.run_grt = runner
    try:
        pygrt.utils.xy2geo(
            "grid_local.nc",
            outgrid="api_grid_geo.nc",
            lat0=35.0,
            lon0=179.5,
        )
        pygrt.utils.geo2xy(
            qfile="local_points.txt",
            outgrid="api_points_local.txt",
            lat0=35.0,
            lon0=179.5,
        )
    finally:
        pygrt.utils.run_grt = original_runner

    for transform in (pygrt.utils.xy2geo, pygrt.utils.geo2xy):
        try:
            transform("input.nc", "output.nc", lat0=35.0, lon0=179.5)
        except TypeError:
            pass
        else:
            raise AssertionError("Only ingrid may be passed positionally.")

    assert runner.commands == [
        ["xy2geo", "-Ggrid_local.nc", "-Oapi_grid_geo.nc", "-C35/179.5"],
        ["geo2xy", "-Qlocal_points.txt", "-Oapi_points_local.txt", "-C35/179.5"],
    ]


assert_text_conversion()
