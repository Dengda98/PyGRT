from pathlib import Path
from tempfile import TemporaryDirectory

import pygrt


MODEL_PARAMS = (8.0, 4.62, 3.3)
COMMON_PARAMS = {
    "modelparams": MODEL_PARAMS,
    "dist": 10.0,
    "nt": 16,
    "dt": 0.01,
    "azimuth": 30.0,
    "scale": 1e20,
    "print_log": False,
}


with TemporaryDirectory(prefix="pygrt-lamb-") as temporary:
    output_root = Path(temporary)

    try:
        pygrt.utils.lamb(
            **COMMON_PARAMS,
            depsrc=5.0,
            deprcv=0.0,
            output_path=output_root / "empty_phases",
            force=(0.5, -1.0, 2.0),
            phases=[],
        )
    except ValueError:
        pass
    else:
        raise ValueError("lamb should reject an empty phase list.")

    pygrt.utils.lamb(
        **COMMON_PARAMS,
        depsrc=0.0,
        deprcv=0.0,
        output_path=output_root / "surface",
        force=(0.5, -1.0, 2.0),
    )
    pygrt.utils.lamb(
        **COMMON_PARAMS,
        depsrc=5.0,
        deprcv=0.0,
        output_path=output_root / "lamb2",
        moment_tensor=(1.0, -2.0, 3.0, 0.5, 1.2, -0.7),
    )
    pygrt.utils.lamb(
        **COMMON_PARAMS,
        depsrc=5.0,
        deprcv=0.0,
        output_path=output_root / "lamb2_phases",
        moment_tensor=(1.0, -2.0, 3.0, 0.5, 1.2, -0.7),
        phases=["P", "S", "SP"],
    )
    pygrt.utils.lamb(
        **COMMON_PARAMS,
        depsrc=5.0,
        deprcv=1.0,
        output_path=output_root / "lamb3",
        strike=100.0,
        dip=30.0,
        rake=70.0,
    )
