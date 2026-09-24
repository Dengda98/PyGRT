import subprocess

import numpy as np
import pygrt


try:
    pygrt.utils.lamb3(nu=0.5, tbar=np.asarray([0.0]), R=10.0, depsrc=1.0, deprcv=1.0, azimuth=0.0)
except ValueError:
    pass
else:
    raise ValueError("lamb3 should reject a Poisson ratio outside (0, 0.5).")


try:
    pygrt.utils.lamb3(nu=0.25, tbar=np.asarray([0.0]), R=0.0, depsrc=1.0, deprcv=1.0, azimuth=0.0)
except ValueError:
    pass
else:
    raise ValueError("lamb3 should require a positive horizontal distance.")


for depsrc, deprcv in ((0.0, 1.0), (1.0, 0.0), (0.0, 0.0)):
    try:
        pygrt.utils.lamb3(
            nu=0.25, tbar=np.asarray([0.0]), R=10.0,
            depsrc=depsrc, deprcv=deprcv, azimuth=0.0,
        )
    except ValueError:
        pass
    else:
        raise ValueError("lamb3 should require both source and receiver below the free surface.")


invalid_lamb3_inputs = (
    ("an empty phase list", dict(nu=0.25, tbar=np.asarray([0.0]), R=10.0, depsrc=2.0, deprcv=1.0, azimuth=0.0, phases=[])),
    ("an empty time series", dict(nu=0.25, tbar=np.asarray([]), R=10.0, depsrc=2.0, deprcv=1.0, azimuth=0.0)),
    ("a multidimensional time series", dict(nu=0.25, tbar=np.asarray([[0.0]]), R=10.0, depsrc=2.0, deprcv=1.0, azimuth=0.0)),
    ("a non-finite time series", dict(nu=0.25, tbar=np.asarray([0.0, np.inf]), R=10.0, depsrc=2.0, deprcv=1.0, azimuth=0.0)),
    ("a non-increasing time series", dict(nu=0.25, tbar=np.asarray([0.0, 0.0]), R=10.0, depsrc=2.0, deprcv=1.0, azimuth=0.0)),
    ("a non-finite Poisson ratio", dict(nu=np.nan, tbar=np.asarray([0.0]), R=10.0, depsrc=2.0, deprcv=1.0, azimuth=0.0)),
    ("a zero horizontal distance", dict(nu=0.25, tbar=np.asarray([0.0]), R=0.0, depsrc=2.0, deprcv=1.0, azimuth=0.0)),
    ("a non-finite horizontal distance", dict(nu=0.25, tbar=np.asarray([0.0]), R=np.inf, depsrc=2.0, deprcv=1.0, azimuth=0.0)),
    ("a non-finite source depth", dict(nu=0.25, tbar=np.asarray([0.0]), R=10.0, depsrc=np.nan, deprcv=1.0, azimuth=0.0)),
    ("a non-finite receiver depth", dict(nu=0.25, tbar=np.asarray([0.0]), R=10.0, depsrc=2.0, deprcv=np.nan, azimuth=0.0)),
    ("a non-finite azimuth", dict(nu=0.25, tbar=np.asarray([0.0]), R=10.0, depsrc=2.0, deprcv=1.0, azimuth=np.nan)),
)
for name, kwargs in invalid_lamb3_inputs:
    try:
        pygrt.utils.lamb3(**kwargs)
    except ValueError:
        pass
    else:
        raise ValueError(f"lamb3 should reject {name}.")


R = 10.0
DEPSRC = 2.0
DEPRCV = 1.0
AZIMUTH = 30.0
# 固定 S 波速度只用于把物理时间换算为无量纲时间
BETA = 1.0
TS = np.arange(0.0, 2.0 + 1e-8, 1e-2)


G, dG_source, dG_receiver, dG_mixed = pygrt.utils.lamb3(
    nu=0.25, tbar=TS, R=R, depsrc=DEPSRC, deprcv=DEPRCV, azimuth=AZIMUTH
)

if G.shape != (len(TS), 3, 3):
    raise ValueError(f"Unexpected lamb3 Green-function shape: {G.shape}")
if dG_source.shape != (len(TS), 3, 3, 3):
    raise ValueError(f"Unexpected lamb3 source-derivative shape: {dG_source.shape}")
if dG_receiver.shape != (len(TS), 3, 3, 3):
    raise ValueError(f"Unexpected lamb3 receiver-derivative shape: {dG_receiver.shape}")
if dG_mixed.shape != (len(TS), 3, 3, 3, 3):
    raise ValueError(f"Unexpected lamb3 mixed-derivative shape: {dG_mixed.shape}")
if not np.isfinite(G).all() or not np.isfinite(dG_source).all() or not np.isfinite(dG_receiver).all() or not np.isfinite(dG_mixed).all():
    raise ValueError("lamb3 returned a non-finite value.")
if not np.allclose(dG_receiver[:, :2], -dG_source[:, :2]):
    raise ValueError("lamb3 horizontal receiver/source derivatives violate translation invariance.")

# 七个震相项彼此独立
phase_results = [
    pygrt.utils.lamb3(
        nu=0.25, tbar=TS, R=R, depsrc=DEPSRC, deprcv=DEPRCV, azimuth=AZIMUTH, phases=[phase]
    )
    for phase in ("P", "S", "PP", "SS", "PS", "SP", "sPs")
]
for index, name in enumerate(("G", "source", "receiver", "mixed")):
    if not np.allclose(
        sum(result[index] for result in phase_results),
        (G, dG_source, dG_receiver, dG_mixed)[index], rtol=1e-9, atol=1e-9,
    ):
        raise ValueError(f"The lamb3 phase-separated {name} results do not add up to the full solution.")

if not np.allclose(
    phase_results[3][0] + phase_results[6][0],
    pygrt.utils.lamb3(
        nu=0.25, tbar=TS, R=R, depsrc=DEPSRC, deprcv=DEPRCV, azimuth=AZIMUTH,
        phases=["SS", "sPs"],
    )[0],
    rtol=1e-9,
    atol=1e-9,
):
    raise ValueError("The lamb3 SS and sPs selections are not independently composable.")

empty = pygrt.utils.lamb3(
    nu=0.25, tbar=TS, R=R, depsrc=DEPSRC, deprcv=DEPRCV, azimuth=AZIMUTH, phases=["BAD"]
)
if any(np.any(values) for values in empty):
    raise ValueError("The lamb3 result should be zero when no valid phase remains.")


def _check_lamb3_right_limit(boundary):
    dt = 1e-2
    epsilon = min(1e-8, dt * 1e-5)
    exact, _, _, _ = pygrt.utils.lamb3(
        nu=0.25, tbar=np.asarray([boundary, boundary + dt]), R=R,
        depsrc=DEPSRC, deprcv=DEPRCV, azimuth=AZIMUTH,
    )
    right, _, _, _ = pygrt.utils.lamb3(
        nu=0.25, tbar=np.asarray([boundary + epsilon]), R=R,
        depsrc=DEPSRC, deprcv=DEPRCV, azimuth=AZIMUTH,
    )
    if not np.allclose(exact[0], right[0], rtol=1e-10, atol=1e-12):
        raise ValueError(f"lamb3 does not use the right-hand value at tbar={boundary:g}.")


main_p_arrival = np.sqrt(0.5 * (1.0 - 2.0 * 0.25) / (1.0 - 0.25))
_check_lamb3_right_limit(main_p_arrival)
_check_lamb3_right_limit(1.0)


reciprocal_closed, reciprocal_closed_source, _, reciprocal_closed_mixed = pygrt.utils.lamb3(
    nu=0.25, tbar=TS, R=R, depsrc=DEPRCV, deprcv=DEPSRC, azimuth=AZIMUTH + 180.0
)
expected_closed_vertical = np.transpose(reciprocal_closed_source[:, 2], (0, 2, 1))
if not np.allclose(dG_receiver[:, 2], expected_closed_vertical, rtol=2e-6, atol=1e-5):
    raise ValueError("lamb3 receiver vertical derivatives violate reciprocity.")
if not np.allclose(dG_mixed, np.transpose(reciprocal_closed_mixed, (0, 2, 1, 4, 3)), rtol=2e-6, atol=1e-5):
    raise ValueError("lamb3 mixed derivatives violate reciprocity.")

for depsrc, deprcv, azimuth in (
    (0.8, 1.7, 0.0),
    (1.7, 0.8, 90.0),
    (4.0, 0.7, 210.0),
    (2.0, 2.0, 45.0),
):
    angle_times = np.asarray([0.7, 1.2, 1.8])
    angle_G, angle_source, angle_receiver, angle_mixed = pygrt.utils.lamb3(
        nu=0.25, tbar=angle_times, R=R, depsrc=depsrc, deprcv=deprcv, azimuth=azimuth
    )
    swapped_G, swapped_source, _, swapped_mixed = pygrt.utils.lamb3(
        nu=0.25, tbar=angle_times, R=R, depsrc=deprcv, deprcv=depsrc,
        azimuth=(azimuth + 180.0) % 360.0,
    )
    if not np.isfinite(angle_G).all() or not np.isfinite(angle_source).all() or not np.isfinite(angle_receiver).all():
        raise ValueError(f"lamb3 returned a non-finite value at case {(depsrc, deprcv, azimuth)}.")
    if not np.isfinite(swapped_G).all() or not np.isfinite(swapped_source).all():
        raise ValueError(f"lamb3 reciprocal case returned a non-finite value at azimuth={azimuth:g}.")
    expected_angle_vertical = np.transpose(swapped_source[:, 2], (0, 2, 1))
    if not np.allclose(angle_receiver[:, 2], expected_angle_vertical, rtol=3e-8, atol=3e-8):
        raise ValueError(f"lamb3 PS/SP vertical reciprocity failed at azimuth={azimuth:g}.")
    if not np.allclose(angle_mixed, np.transpose(swapped_mixed, (0, 2, 1, 4, 3)), rtol=3e-8, atol=3e-8):
        raise ValueError(f"lamb3 mixed reciprocity failed at azimuth={azimuth:g}.")


critical_edge = np.array([1.612])
edge_G, edge_dG_source, edge_dG_receiver, _ = pygrt.utils.lamb3(
    nu=0.25, tbar=critical_edge, R=10.0, depsrc=10.0, deprcv=5.0, azimuth=30.0
)
if not np.isfinite(edge_G).all() or not np.isfinite(edge_dG_source).all() or not np.isfinite(edge_dG_receiver).all():
    raise ValueError("lamb3 failed at the subcritical sPs edge case.")

cli_result = subprocess.run(
    [
        "grt",
        "lamb3",
        "-P0.25",
        "-T0/2/1e-2",
        f"-R{R}",
        f"-D{DEPSRC}/{DEPRCV}",
        "-S+slamb3_source+rlamb3_receiver",
        f"-A{AZIMUTH}",
    ],
    check=True,
    capture_output=True,
    text=True,
)
cli = np.loadtxt(cli_result.stdout.splitlines()[1:])
if cli.shape != (len(TS), 10):
    raise ValueError(f"Unexpected lamb3 CLI shape: {cli.shape}")
if not np.allclose(cli[:, 1:10], G.reshape(len(TS), 9), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb3 CLI and Python Green functions differ.")
source_cli = np.loadtxt("lamb3_source")
receiver_cli = np.loadtxt("lamb3_receiver")
if source_cli.shape != (len(TS), 28) or receiver_cli.shape != (len(TS), 28):
    raise ValueError(f"Unexpected lamb3 derivative CLI shapes: {source_cli.shape}, {receiver_cli.shape}")
if not np.allclose(source_cli[:, 1:], dG_source.reshape(len(TS), 27), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb3 source derivative file and Python derivatives differ.")
if not np.allclose(receiver_cli[:, 1:], dG_receiver.reshape(len(TS), 27), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb3 receiver derivative file and Python derivatives differ.")
mixed_cli = np.loadtxt("lamb3_mixed")
if mixed_cli.shape != (len(TS), 82) or not np.isfinite(mixed_cli).all():
    raise ValueError(f"Unexpected lamb3 mixed derivative file shape or values: {mixed_cli.shape}")
if not np.allclose(mixed_cli[:, 1:], dG_mixed.reshape(len(TS), 81), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb3 mixed derivative file and Python derivatives differ.")

phase_cli = np.loadtxt("lamb3_phases")[:, 1:].reshape(len(TS), 3, 3)
phase_python = pygrt.utils.lamb3(
    nu=0.25, tbar=TS, R=R, depsrc=DEPSRC, deprcv=DEPRCV, azimuth=AZIMUTH,
    phases=["P", "S", "PP", "SS", "PS", "SP", "sPs"],
)[0]
if not np.allclose(phase_cli, phase_python, rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb3 phase-list CLI and Python results differ.")


def _lamb3_geometry(source, receiver):
    horizontal = receiver[:2] - source[:2]
    horizontal_distance = np.hypot(horizontal[0], horizontal[1])
    azimuth = np.degrees(np.arctan2(horizontal[1], horizontal[0])) % 360.0
    return horizontal_distance, source[2], receiver[2], azimuth


def _lamb3_green_at_physical_time(physical_time, source, receiver):
    distance = np.linalg.norm(receiver - source)
    horizontal_distance, depsrc, deprcv, azimuth = _lamb3_geometry(source, receiver)
    return pygrt.utils.lamb3(
        nu=0.25,
        tbar=np.asarray([physical_time * BETA / distance]),
        R=horizontal_distance,
        depsrc=depsrc,
        deprcv=deprcv,
        azimuth=azimuth,
    )[0][0]


def _check_lamb3_mixed_geometry_derivatives():
    source = np.array([0.0, 0.0, 2.0])
    receiver = np.array([
        R * np.cos(np.deg2rad(AZIMUTH)),
        R * np.sin(np.deg2rad(AZIMUTH)),
        1.0,
    ])
    distance = np.linalg.norm(receiver - source)
    center = 1.8
    physical_time = center * distance / BETA
    _, _, _, expected_series = pygrt.utils.lamb3(
        nu=0.25,
        tbar=np.asarray([center - 2e-3, center, center + 2e-3]),
        R=R,
        depsrc=source[2],
        deprcv=receiver[2],
        azimuth=AZIMUTH,
    )
    expected = expected_series[1]
    finite_difference = np.empty_like(expected)
    step = 1e-2
    for receiver_direction in range(3):
        for source_direction in range(3):
            source_plus = source.copy()
            source_minus = source.copy()
            receiver_plus = receiver.copy()
            receiver_minus = receiver.copy()
            source_plus[source_direction] += step
            source_minus[source_direction] -= step
            receiver_plus[receiver_direction] += step
            receiver_minus[receiver_direction] -= step
            values = (
                _lamb3_green_at_physical_time(physical_time, source_plus, receiver_plus)
                / np.linalg.norm(receiver_plus - source_plus),
                _lamb3_green_at_physical_time(physical_time, source_plus, receiver_minus)
                / np.linalg.norm(receiver_minus - source_plus),
                _lamb3_green_at_physical_time(physical_time, source_minus, receiver_plus)
                / np.linalg.norm(receiver_plus - source_minus),
                _lamb3_green_at_physical_time(physical_time, source_minus, receiver_minus)
                / np.linalg.norm(receiver_minus - source_minus),
            )
            finite_difference[receiver_direction, source_direction] = (
                values[0] - values[1] - values[2] + values[3]
            ) * distance**3 / (4.0 * step**2)
    if not np.allclose(finite_difference, expected, rtol=5e-3, atol=2e-2):
        raise ValueError("lamb3 mixed derivatives fail the geometry finite-difference test.")


def _check_lamb3_geometry_derivatives(horizontal_distance, depsrc, deprcv):
    source = np.array([0.0, 0.0, depsrc])
    receiver = np.array([
        horizontal_distance * np.cos(np.deg2rad(AZIMUTH)),
        horizontal_distance * np.sin(np.deg2rad(AZIMUTH)),
        deprcv,
    ])
    distance = np.linalg.norm(receiver - source)
    times = (0.70, 0.85, 1.20, 1.50, 1.80)
    for tbar in times:
        _, expected_source, expected_receiver, _ = pygrt.utils.lamb3(
            nu=0.25,
            tbar=np.asarray([tbar - 1e-3, tbar, tbar + 1e-3]),
            R=horizontal_distance,
            depsrc=depsrc,
            deprcv=deprcv,
            azimuth=AZIMUTH,
        )
        expected_source = expected_source[1]
        expected_receiver = expected_receiver[1]
        physical_time = tbar * distance / BETA
        for kind, expected in (("source", expected_source), ("receiver", expected_receiver)):
            for coordinate in range(3):
                source_plus = source.copy()
                source_minus = source.copy()
                receiver_plus = receiver.copy()
                receiver_minus = receiver.copy()
                if kind == "source":
                    source_plus[coordinate] += 1e-2
                    source_minus[coordinate] -= 1e-2
                else:
                    receiver_plus[coordinate] += 1e-2
                    receiver_minus[coordinate] -= 1e-2
                distance_plus = np.linalg.norm(receiver_plus - source_plus)
                distance_minus = np.linalg.norm(receiver_minus - source_minus)
                finite_difference = (
                    _lamb3_green_at_physical_time(physical_time, source_plus, receiver_plus) / distance_plus
                    - _lamb3_green_at_physical_time(physical_time, source_minus, receiver_minus) / distance_minus
                ) * distance**2 / (2e-2)
                if not np.allclose(finite_difference, expected[coordinate], rtol=5e-3, atol=2e-3):
                    raise ValueError(
                        f"lamb3 {kind} derivative failed at case "
                        f"({horizontal_distance:g}, {depsrc:g}, {deprcv:g}), "
                        f"tbar={tbar:g}, coordinate={coordinate + 1}"
                    )


for depsrc, deprcv in ((0.5, 0.1), (2.0, 1.0), (10.0, 5.0)):
    _check_lamb3_geometry_derivatives(10.0, depsrc, deprcv)
_check_lamb3_mixed_geometry_derivatives()

shallow_ts = np.arange(1.49, 1.53 + 1e-8, 2e-3)
_, shallow_source, shallow_receiver, _ = pygrt.utils.lamb3(
    nu=0.25, tbar=shallow_ts, R=10.0, depsrc=0.5, deprcv=0.1, azimuth=30.0
)
for name, values in (("source", shallow_source), ("receiver", shallow_receiver)):
    jump = np.max(np.abs(np.diff(values[:, 2, 2, 2])))
    if jump > 1e-2:
        raise ValueError(f"lamb3 {name} G33,3 has a nonphysical jump near tbar=1.5: {jump:g}")

print("lamb3 interface, geometry-derivative, mixed, and continuity tests passed.")
