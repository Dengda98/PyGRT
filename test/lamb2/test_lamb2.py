import subprocess

import numpy as np
import pygrt


def expect_value_error(desc, **kwargs):
    try:
        pygrt.utils.lamb2(**kwargs)
    except ValueError:
        return
    raise ValueError(f"lamb2 should reject {desc}.")


expect_value_error("a Poisson ratio outside (0, 0.5)", nu=0.0, tbar=np.asarray([0.0]), R=10.0, depsrc=5.0, azimuth=0.0)
expect_value_error("a non-positive depsrc", nu=0.25, tbar=np.asarray([0.0]), R=10.0, depsrc=0.0, azimuth=0.0)
expect_value_error("a non-positive deprcv", nu=0.25, tbar=np.asarray([0.0]), R=10.0, deprcv=0.0, azimuth=0.0)
expect_value_error("a negative depsrc", nu=0.25, tbar=np.asarray([0.0]), R=10.0, depsrc=-1.0, azimuth=0.0)
expect_value_error("a negative deprcv", nu=0.25, tbar=np.asarray([0.0]), R=10.0, deprcv=-1.0, azimuth=0.0)
expect_value_error("both depsrc and deprcv", nu=0.25, tbar=np.asarray([0.0]), R=10.0, depsrc=5.0, deprcv=3.0, azimuth=0.0)
expect_value_error("neither depsrc nor deprcv", nu=0.25, tbar=np.asarray([0.0]), R=10.0, azimuth=0.0)
expect_value_error("a non-positive horizontal distance", nu=0.25, tbar=np.asarray([0.0]), R=0.0, depsrc=5.0, azimuth=0.0)

invalid_lamb2_inputs = (
    ("an empty time series", dict(nu=0.25, tbar=np.asarray([]), R=10.0, depsrc=5.0, azimuth=0.0)),
    ("a multidimensional time series", dict(nu=0.25, tbar=np.asarray([[0.0]]), R=10.0, depsrc=5.0, azimuth=0.0)),
    ("a non-finite time series", dict(nu=0.25, tbar=np.asarray([0.0, np.nan]), R=10.0, depsrc=5.0, azimuth=0.0)),
    ("a non-increasing time series", dict(nu=0.25, tbar=np.asarray([0.0, 0.0]), R=10.0, depsrc=5.0, azimuth=0.0)),
    ("a non-finite Poisson ratio", dict(nu=np.nan, tbar=np.asarray([0.0]), R=10.0, depsrc=5.0, azimuth=0.0)),
    ("a negative horizontal distance", dict(nu=0.25, tbar=np.asarray([0.0]), R=-1.0, depsrc=5.0, azimuth=0.0)),
    ("a zero horizontal distance", dict(nu=0.25, tbar=np.asarray([0.0]), R=0.0, depsrc=5.0, azimuth=0.0)),
    ("a non-finite horizontal distance", dict(nu=0.25, tbar=np.asarray([0.0]), R=np.inf, depsrc=5.0, azimuth=0.0)),
    ("a non-finite depsrc", dict(nu=0.25, tbar=np.asarray([0.0]), R=10.0, depsrc=np.nan, azimuth=0.0)),
    ("a non-finite deprcv", dict(nu=0.25, tbar=np.asarray([0.0]), R=10.0, deprcv=np.nan, azimuth=0.0)),
    ("a non-finite azimuth", dict(nu=0.25, tbar=np.asarray([0.0]), R=10.0, depsrc=5.0, azimuth=np.nan)),
)
for name, kwargs in invalid_lamb2_inputs:
    expect_value_error(name, **kwargs)


R = 10.0
SOURCE_DEPTH = 5.0
AZIMUTH = 30.0
ts = np.arange(0.0, 2.0 + 1e-8, 1e-2)
G, dG_source, dG_receiver = pygrt.utils.lamb2(nu=0.25, tbar=ts, R=R, depsrc=SOURCE_DEPTH, azimuth=AZIMUTH)

if not np.allclose(dG_receiver[:, :2], -dG_source[:, :2]):
    raise ValueError("Horizontal receiver and source derivatives violate translation invariance.")
if not np.allclose(dG_receiver[:, 2, 0, :], dG_source[:, 0, 2, :]):
    raise ValueError("The vertical receiver derivative first row is inconsistent.")
if not np.allclose(dG_receiver[:, 2, 1, 0], dG_source[:, 0, 2, 1]):
    raise ValueError("The vertical receiver derivative (2, 1) component is inconsistent.")
if not np.allclose(dG_receiver[:, 2, 1, 1:], dG_source[:, 1, 2, 1:]):
    raise ValueError("The vertical receiver derivative second row is inconsistent.")

for azimuth in (0.0, 90.0, 180.0, 270.0):
    angle_G, angle_source, angle_receiver = pygrt.utils.lamb2(
        nu=0.25, tbar=np.asarray([0.0, 0.65, 1.05, 1.8]), R=R, depsrc=SOURCE_DEPTH, azimuth=azimuth
    )
    if not np.isfinite(angle_G).all() or not np.isfinite(angle_source).all() or not np.isfinite(angle_receiver).all():
        raise ValueError(f"lamb2 returned a non-finite value at azimuth={azimuth:g}.")
    if not np.allclose(angle_receiver[:, :2], -angle_source[:, :2]):
        raise ValueError(f"lamb2 horizontal derivative relation failed at azimuth={azimuth:g}.")


def _check_reciprocity(depth, azimuth):
    buried_G, buried_source, buried_receiver = pygrt.utils.lamb2(
        nu=0.25, tbar=ts, R=R, depsrc=depth, azimuth=azimuth
    )
    surface_G, surface_source, surface_receiver = pygrt.utils.lamb2(
        nu=0.25, tbar=ts, R=R, deprcv=depth, azimuth=(azimuth + 180.0) % 360.0
    )
    if not np.allclose(surface_G, np.swapaxes(buried_G, -2, -1), rtol=2e-6, atol=1e-5):
        raise ValueError(f"lamb2 Green-function reciprocity failed at azimuth={azimuth:g}.")
    if not np.allclose(surface_source, np.swapaxes(buried_receiver, -2, -1), rtol=2e-6, atol=1e-5):
        raise ValueError(f"lamb2 source-derivative reciprocity failed at azimuth={azimuth:g}.")
    if not np.allclose(surface_receiver, np.swapaxes(buried_source, -2, -1), rtol=2e-6, atol=1e-5):
        raise ValueError(f"lamb2 receiver-derivative reciprocity failed at azimuth={azimuth:g}.")
    if not np.allclose(surface_receiver[:, :2], -surface_source[:, :2]):
        raise ValueError(f"lamb2 surface-source horizontal derivatives violate translation invariance at azimuth={azimuth:g}.")


for azimuth in (0.0, 30.0, 90.0, 180.0, 210.0, 360.0):
    _check_reciprocity(SOURCE_DEPTH, azimuth)
for depth in (0.5, 2.0, 10.0):
    _check_reciprocity(depth, AZIMUTH)


G_surface, dG_surface_source, dG_surface_receiver = pygrt.utils.lamb2(
    nu=0.25, tbar=ts, R=R, deprcv=SOURCE_DEPTH, azimuth=AZIMUTH
)
if not np.isfinite(G_surface).all() or not np.isfinite(dG_surface_source).all() or not np.isfinite(dG_surface_receiver).all():
    raise ValueError("lamb2 returned a non-finite value for the surface-source case.")
if not np.allclose(dG_surface_receiver[:, :2], -dG_surface_source[:, :2]):
    raise ValueError("Surface-source horizontal receiver and source derivatives violate translation invariance.")

cli_result = subprocess.run(
    [
        "grt",
        "lamb2",
        "-P0.25",
        "-T0/2/1e-2",
        f"-R{R}",
        f"-Ds{SOURCE_DEPTH}",
        "-S+slamb2_source+rlamb2_receiver",
        f"-A{AZIMUTH}",
    ],
    check=True,
    capture_output=True,
    text=True,
)
cli = np.loadtxt(cli_result.stdout.splitlines()[1:])

if cli.shape != (len(ts), 10):
    raise ValueError(f"Unexpected lamb2 CLI shape: {cli.shape}")
if not np.allclose(cli[:, 1:10], G.reshape(len(ts), 9), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb2 CLI and Python Green functions differ.")
source_cli = np.loadtxt("lamb2_source")
receiver_cli = np.loadtxt("lamb2_receiver")
if source_cli.shape != (len(ts), 28) or receiver_cli.shape != (len(ts), 28):
    raise ValueError(f"Unexpected lamb2 derivative CLI shapes: {source_cli.shape}, {receiver_cli.shape}")
if not np.allclose(source_cli[:, 1:], dG_source.reshape(len(ts), 27), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb2 source derivative file and Python derivatives differ.")
if not np.allclose(receiver_cli[:, 1:], dG_receiver.reshape(len(ts), 27), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb2 receiver derivative file and Python derivatives differ.")

surface_cli_result = subprocess.run(
    [
        "grt",
        "lamb2",
        "-P0.25",
        "-T0/2/1e-2",
        f"-R{R}",
        f"-Dr{SOURCE_DEPTH}",
        "-S+slamb2_surface_source+rlamb2_surface_receiver",
        f"-A{AZIMUTH}",
    ],
    check=True,
    capture_output=True,
    text=True,
)
surface_cli = np.loadtxt(surface_cli_result.stdout.splitlines()[1:])
if surface_cli.shape != (len(ts), 10):
    raise ValueError(f"Unexpected lamb2 surface-source CLI shape: {surface_cli.shape}")
if not np.allclose(surface_cli[:, 1:10], G_surface.reshape(len(ts), 9), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb2 surface-source CLI and Python Green functions differ.")
surface_source_cli = np.loadtxt("lamb2_surface_source")
surface_receiver_cli = np.loadtxt("lamb2_surface_receiver")
if surface_source_cli.shape != (len(ts), 28) or surface_receiver_cli.shape != (len(ts), 28):
    raise ValueError(
        f"Unexpected lamb2 surface-source derivative CLI shapes: {surface_source_cli.shape}, {surface_receiver_cli.shape}"
    )
if not np.allclose(surface_source_cli[:, 1:], dG_surface_source.reshape(len(ts), 27), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb2 surface-source derivative file and Python derivatives differ.")
if not np.allclose(surface_receiver_cli[:, 1:], dG_surface_receiver.reshape(len(ts), 27), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb2 surface-receiver derivative file and Python derivatives differ.")


def _check_lamb2_right_limit(boundary, **depth_kw):
    dt = 1e-2
    epsilon = min(1e-8, dt * 1e-5)
    exact, _, _ = pygrt.utils.lamb2(
        nu=0.25, tbar=np.asarray([boundary, boundary + dt]), R=R, azimuth=AZIMUTH, **depth_kw
    )
    right, _, _ = pygrt.utils.lamb2(
        nu=0.25, tbar=np.asarray([boundary + epsilon]), R=R, azimuth=AZIMUTH, **depth_kw
    )
    if not np.allclose(exact[0], right[0], rtol=1e-10, atol=1e-12):
        raise ValueError(f"lamb2 does not use the right-hand value at tbar={boundary:g}.")


main_p_arrival = np.sqrt(0.5 * (1.0 - 2.0 * 0.25) / (1.0 - 0.25))
for depth_kw in ({"depsrc": SOURCE_DEPTH}, {"deprcv": SOURCE_DEPTH}):
    _check_lamb2_right_limit(main_p_arrival, **depth_kw)
    _check_lamb2_right_limit(1.0, **depth_kw)
    main_theta = np.arctan2(R, SOURCE_DEPTH)
    main_t_sp = np.cos(main_theta - np.arcsin(main_p_arrival))
    _check_lamb2_right_limit(main_t_sp, **depth_kw)


def _lamb2_g_over_r(nu, tbar, source, reference_radius):
    """Return ``Gbar / r`` at fixed physical times for a buried source and a surface receiver."""
    source_to_receiver = -np.asarray(source, dtype=float)
    radius = np.linalg.norm(source_to_receiver)
    ray = source_to_receiver / radius
    horizontal_distance = np.hypot(source_to_receiver[0], source_to_receiver[1])
    source_depth = source[2]
    azimuth = np.rad2deg(np.arctan2(ray[1], ray[0])) % 360.0
    G, _, _ = pygrt.utils.lamb2(
        nu=nu, tbar=tbar * reference_radius / radius, R=horizontal_distance,
        depsrc=source_depth, azimuth=azimuth,
    )
    return G / radius


fd_ts = np.arange(0.05, 1.95 + 0.5e-3, 1e-3)
theta_rad = np.deg2rad(60.0)
azimuth_rad = np.deg2rad(30.0)
source = np.array([
    -np.sin(theta_rad) * np.cos(azimuth_rad),
    -np.sin(theta_rad) * np.sin(azimuth_rad),
    np.cos(theta_rad),
])
reference_radius = np.linalg.norm(source)
_, dG_source_fd, _ = pygrt.utils.lamb2(
    nu=0.25, tbar=fd_ts, R=np.sin(theta_rad), depsrc=np.cos(theta_rad), azimuth=30.0
)
p_arrival = np.sqrt(0.5 * (1.0 - 2.0 * 0.25) / (1.0 - 0.25))
theta_c = np.arcsin(p_arrival)
arrival_times = (p_arrival, 1.0, np.cos(theta_rad - theta_c))
smooth = np.ones(fd_ts.shape, dtype=bool)
for arrival_time in arrival_times:
    smooth &= np.abs(fd_ts - arrival_time) > 0.03

target = dG_source_fd.transpose(1, 0, 2, 3)
relative_errors = []
for step in (1e-2, 3e-3, 1e-3):
    finite_difference = np.empty_like(target)
    for coordinate in range(3):
        source_plus = source.copy()
        source_minus = source.copy()
        source_plus[coordinate] += step
        source_minus[coordinate] -= step
        finite_difference[coordinate] = (
            _lamb2_g_over_r(0.25, fd_ts, source_plus, reference_radius)
            - _lamb2_g_over_r(0.25, fd_ts, source_minus, reference_radius)
        ) * reference_radius**2 / (2.0 * step)

    error = np.abs(finite_difference[:, smooth] - target[:, smooth])
    scale = np.max(np.abs(target[:, smooth]))
    relative_error = np.max(error) / scale
    relative_errors.append(relative_error)
    print(
        f"lamb2 source finite-difference result for step={step:g}: "
        f"max error={np.max(error):.3e}, relative error={relative_error:.3e}."
    )

if not relative_errors[0] > relative_errors[1] > relative_errors[2]:
    raise ValueError(f"The lamb2 source finite-difference error does not converge: {relative_errors}.")
if not np.allclose(finite_difference[:, smooth], target[:, smooth], rtol=5e-4, atol=1e-3):
    raise ValueError(
        f"The lamb2 source derivatives fail the source-position finite-difference test: "
        f"relative error={relative_errors[-1]:.3e}."
    )
print(f"lamb2 source finite-difference test passed: relative error={relative_errors[-1]:.3e}.")


def _lamb2_geometry(source, receiver):
    horizontal = receiver[:2] - source[:2]
    horizontal_distance = np.hypot(horizontal[0], horizontal[1])
    azimuth = np.degrees(np.arctan2(horizontal[1], horizontal[0])) % 360.0
    return horizontal_distance, source[2], receiver[2], azimuth


def _lamb2_green_at_physical_time(physical_time, source, receiver):
    distance = np.linalg.norm(receiver - source)
    horizontal_distance, source_depth, receiver_depth, azimuth = _lamb2_geometry(source, receiver)
    if source_depth > 0.0 and receiver_depth == 0.0:
        depth_kw = {"depsrc": source_depth}
    elif source_depth == 0.0 and receiver_depth > 0.0:
        depth_kw = {"deprcv": receiver_depth}
    else:
        raise ValueError("lamb2 finite-difference geometry must keep exactly one point underground.")
    return pygrt.utils.lamb2(
        nu=0.25,
        tbar=np.asarray([physical_time * 1.0 / distance]),
        R=horizontal_distance,
        azimuth=azimuth,
        **depth_kw,
    )[0][0]


def _check_surface_source_receiver_derivatives():
    source = np.array([0.0, 0.0, 0.0])
    receiver = np.array([
        R * np.cos(np.deg2rad(AZIMUTH)),
        R * np.sin(np.deg2rad(AZIMUTH)),
        SOURCE_DEPTH,
    ])
    distance = np.linalg.norm(receiver - source)
    times = (0.70, 1.20, 1.80)
    for tbar in times:
        _, _, expected_receiver = pygrt.utils.lamb2(
            nu=0.25,
            tbar=np.asarray([tbar - 1e-3, tbar, tbar + 1e-3]),
            R=R,
            deprcv=SOURCE_DEPTH,
            azimuth=AZIMUTH,
        )
        expected_receiver = expected_receiver[1]
        physical_time = tbar * distance
        for coordinate in range(3):
            receiver_plus = receiver.copy()
            receiver_minus = receiver.copy()
            receiver_plus[coordinate] += 1e-2
            receiver_minus[coordinate] -= 1e-2
            if receiver_plus[2] <= 0.0 or receiver_minus[2] <= 0.0:
                raise ValueError("Surface-source receiver finite difference left the halfspace.")
            distance_plus = np.linalg.norm(receiver_plus - source)
            distance_minus = np.linalg.norm(receiver_minus - source)
            finite_difference = (
                _lamb2_green_at_physical_time(physical_time, source, receiver_plus) / distance_plus
                - _lamb2_green_at_physical_time(physical_time, source, receiver_minus) / distance_minus
            ) * distance**2 / (2e-2)
            if not np.allclose(finite_difference, expected_receiver[coordinate], rtol=5e-3, atol=2e-3):
                raise ValueError(
                    f"lamb2 surface-source receiver derivative failed at tbar={tbar:g}, "
                    f"coordinate={coordinate + 1}"
                )


_check_surface_source_receiver_derivatives()
print("lamb2 surface-source receiver finite-difference test passed.")

print("lamb2 C/Python interface test passed.")
