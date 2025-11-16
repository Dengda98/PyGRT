import numpy as np
import pygrt

ts = np.arange(0, 2+1e-8, 1e-3)
lamb1 = pygrt.utils.solve_lamb1(0.25, ts, 30)

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