from pathlib import Path

from pygrt.utils import compute_okada


def main():
    Path("cfaults.inp").write_text(
        "  #   X-start    Y-start     X-fin     Y-fin  Kode  shear(m)  reverse(m)  dip  top  bot\n"
        "xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx  xxxxxxxxxx  xxx  xxx  xxx\n"
        "  1     0.0000     0.0000    30.0000     0.0000 100     0.1000      0.2000  70.0 20.0 30.0\n",
        encoding="utf-8",
    )

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
    Path("cfaults.inp").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
