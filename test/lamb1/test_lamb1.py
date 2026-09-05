import numpy as np
import pygrt


def expect_value_error(desc, **kwargs):
    try:
        pygrt.utils.lamb1(**kwargs)
    except ValueError:
        return
    raise ValueError(f"lamb1 should reject {desc}.")


expect_value_error("a Poisson ratio outside (0, 0.5)", nu=0.5, tbar=np.asarray([0.0]), azimuth=0.0)
expect_value_error("an empty time series", nu=0.25, tbar=np.asarray([]), azimuth=0.0)
expect_value_error("a multidimensional time series", nu=0.25, tbar=np.asarray([[0.0]]), azimuth=0.0)
expect_value_error("a non-finite time series", nu=0.25, tbar=np.asarray([0.0, np.inf]), azimuth=0.0)
expect_value_error("a negative time series", nu=0.25, tbar=np.asarray([-1e-8]), azimuth=0.0)
expect_value_error("a non-increasing time series", nu=0.25, tbar=np.asarray([0.0, 0.0]), azimuth=0.0)
expect_value_error("a non-finite Poisson ratio", nu=np.nan, tbar=np.asarray([0.0]), azimuth=0.0)
expect_value_error("a non-finite azimuth", nu=0.25, tbar=np.asarray([0.0]), azimuth=np.nan)
expect_value_error("an azimuth outside [0, 360]", nu=0.25, tbar=np.asarray([0.0]), azimuth=361.0)


ts = np.arange(0, 2+1e-8, 1e-3)
lamb1 = pygrt.utils.lamb1(nu=0.25, tbar=ts, azimuth=30)

# check results
# Since the rounding error of ref_lamb1, the error will not be absolutely zero.
ref_lamb1 = np.loadtxt("ref_lamb1")[:, 1:].reshape(-1, 3, 3)
tol = 0.01
err = np.sum(np.abs(ref_lamb1 - lamb1)) / np.mean(np.abs(ref_lamb1))
print(f"err={err:e}")
if err > tol:
    raise ValueError(f"err({err}) > tol")

# check c results
c_lamb1 = np.loadtxt("lamb1")[:, 1:].reshape(-1, 3, 3)
err = np.sum(np.abs(ref_lamb1 - c_lamb1)) / np.mean(np.abs(ref_lamb1))
print(f"err={err:e}")
if err > tol:
    raise ValueError(f"err({err}) > tol")
