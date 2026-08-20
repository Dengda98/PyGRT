from pathlib import Path

from pygrt.utils import compute_okada


def main():
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
        output_path="okada_python_ff.nc",
        finite_fault="cfaults.inp",
    )

    print("test_okada.py: all checks passed")

    Path("okada_python.nc").unlink(missing_ok=True)
    Path("okada_python_ff.nc").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
