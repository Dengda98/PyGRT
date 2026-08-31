import subprocess

import numpy as np
import pygrt


ts = np.arange(0.0, 2.0 + 1e-8, 1e-2)
G, dG_source, dG_receiver = pygrt.utils.lamb2(0.25, ts, 60.0, 30.0)

if not np.allclose(dG_receiver[:, :2], -dG_source[:, :2]):
    raise ValueError("Horizontal receiver and source derivatives violate translation invariance.")
if not np.allclose(dG_receiver[:, 2, 0, :], dG_source[:, 0, 2, :]):
    raise ValueError("The vertical receiver derivative first row is inconsistent.")
if not np.allclose(dG_receiver[:, 2, 1, 0], dG_source[:, 0, 2, 1]):
    raise ValueError("The vertical receiver derivative (2, 1) component is inconsistent.")
if not np.allclose(dG_receiver[:, 2, 1, 1:], dG_source[:, 1, 2, 1:]):
    raise ValueError("The vertical receiver derivative second row is inconsistent.")

cli_result = subprocess.run(
    ["grt", "lamb2", "-P0.25", "-T0/2/1e-2", "-D60", "-A30"],
    check=True,
    capture_output=True,
    text=True,
)
cli = np.loadtxt(cli_result.stdout.splitlines()[1:])

if cli.shape != (len(ts), 64):
    raise ValueError(f"Unexpected lamb2 CLI shape: {cli.shape}")
if not np.allclose(cli[:, 1:10], G.reshape(len(ts), 9), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb2 CLI and Python Green functions differ.")
if not np.allclose(cli[:, 10:37], dG_source.reshape(len(ts), 27), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb2 source CLI and Python derivatives differ.")
if not np.allclose(cli[:, 37:], dG_receiver.reshape(len(ts), 27), rtol=2e-6, atol=1e-5):
    raise ValueError("The lamb2 receiver CLI and Python derivatives differ.")

print("lamb2 C/Python interface test passed.")
