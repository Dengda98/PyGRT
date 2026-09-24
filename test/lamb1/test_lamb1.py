import numpy as np
import pygrt


nu = 0.25
azimuth = 30.0
tbar = np.arange(0.0, 2.0 + 1e-8, 1e-3)
cbar = 2e-4


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


for invalid_cbar in (0.0, -1e-6):
    try:
        pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth, cbar=invalid_cbar)
    except ValueError:
        continue
    raise RuntimeError(f"lamb1 should reject cbar={invalid_cbar}.")

try:
    pygrt.utils.lamb1(nu=nu, tbar=np.asarray([0.0]), azimuth=azimuth, phases=[])
except ValueError:
    pass
else:
    raise RuntimeError("lamb1 should reject an empty phase list.")


fixed = pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth)
moving = pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth, cbar=cbar)
require(fixed.shape == (len(tbar), 3, 3), "fixed-source lamb1 should return shape (nt, 3, 3)")
require(moving.shape == (len(tbar), 3), "moving-source lamb1 should return shape (nt, 3)")
require(np.all(np.isfinite(fixed)), "fixed-source lamb1 result should be finite")
require(np.all(np.isfinite(moving)), "moving-source lamb1 result should be finite")

fixed_cli = np.loadtxt("lamb1")
moving_cli = np.loadtxt("lamb1_moving")
require(fixed_cli.shape == (len(tbar), 10), "fixed-source CLI should output 10 columns")
require(moving_cli.shape == (len(tbar), 4), "moving-source CLI should output 4 columns")
require(np.allclose(fixed_cli[:, 0], tbar, rtol=0.0, atol=5e-7), "CLI time column should match tbar")
require(np.allclose(moving_cli[:, 0], tbar, rtol=0.0, atol=5e-7), "moving CLI time column should match tbar")
require(np.allclose(fixed_cli[:, 1:].reshape(-1, 3, 3), fixed, rtol=1e-6, atol=1e-6),
        "fixed-source CLI and Python results should agree")
require(np.allclose(moving_cli[:, 1:], moving, rtol=1e-6, atol=1e-6),
        "moving-source CLI and Python results should agree")

phase_sum = sum(
    pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth, phases=[phase])
    for phase in ("P", "S", "R")
)
require(np.allclose(phase_sum, fixed, rtol=1e-10, atol=1e-10),
        "fixed-source phase-separated results should add up to the full solution")

empty = pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth, phases=["BAD"])
require(np.all(empty == 0.0), "lamb1 should return zeros when no valid phase remains")

phase_cli = np.loadtxt("lamb1_PS")[:, 1:].reshape(-1, 3, 3)
phase_python = pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth, phases=["P", "S"])
require(np.allclose(phase_cli, phase_python, rtol=2e-6, atol=1e-5),
        "lamb1 phase-list CLI and Python results should agree")

rayleigh_cli = np.loadtxt("lamb1_R")[:, 1:].reshape(-1, 3, 3)
rayleigh_python = pygrt.utils.lamb1(nu=nu, tbar=tbar, azimuth=azimuth, phases="R")
require(np.allclose(rayleigh_cli, rayleigh_python, rtol=2e-6, atol=1e-5),
        "lamb1 Rayleigh phase CLI and Python results should agree")

print("lamb1 tests passed")
