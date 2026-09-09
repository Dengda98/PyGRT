import shutil
from pathlib import Path

import numpy as np
from obspy import read

import pygrt


MODEL = "../milrow"


def compare_sac_dirs(c_dir: Path, py_dir: Path, expected_names: set[str]) -> None:
    c_names = {path.name for path in c_dir.glob("*.sac")}
    py_names = {path.name for path in py_dir.glob("*.sac")}
    assert c_names == expected_names
    assert py_names == expected_names

    for name in sorted(expected_names):
        c_trace = read(str(c_dir / name))[0]
        py_trace = read(str(py_dir / name))[0]
        assert c_trace.stats.npts == py_trace.stats.npts
        assert c_trace.stats.delta == py_trace.stats.delta
        np.testing.assert_array_equal(c_trace.data, py_trace.data)


def expect_value_error(callable_obj, message: str) -> None:
    try:
        callable_obj()
    except ValueError:
        return
    raise AssertionError(message)


model = pygrt.PyModel1D(modelpath=MODEL)

# 与 test_rftn.sh 中的 C 命令对应的 P 波参数
model.rftn(
    wtype="P",
    rayp=0.12,
    nt=64,
    dt=0.05,
    output_path="PY_P",
    write_components=True,
    print_log=False,
)
compare_sac_dirs(
    Path("C_P"),
    Path("PY_P"),
    {"P_rftn.sac", "P_Z.sac", "P_R.sac"},
)

# 与 test_rftn.sh 中的 C 命令对应的 SV 波参数
model.rftn(
    wtype="S",
    inca=20.0,
    idx=2,
    nt=64,
    dt=0.05,
    output_path="PY_S",
    zeta=0.9,
    upsampling_n=2,
    keepAllFreq=True,
    skipImagComps=True,
    alp=0.8,
    delay=0.5,
    write_components=True,
    print_log=False,
)
compare_sac_dirs(
    Path("C_S"),
    Path("PY_S"),
    {"S_rftn.sac", "S_Z.sac", "S_R.sac"},
)

# Python 端输入检查
expect_value_error(
    lambda: model.rftn(wtype="P", nt=8, dt=0.1, output_path="bad", rayp=0.1, inca=20.0),
    "rayp and inca should be mutually exclusive",
)
expect_value_error(
    lambda: model.rftn(wtype="P", nt=8, dt=0.1, output_path="bad"),
    "one of rayp and inca should be required",
)
expect_value_error(
    lambda: model.rftn(wtype="P", nt=8, dt=0.1, output_path="bad", rayp=0.1, idx=1),
    "idx should require inca",
)
expect_value_error(
    lambda: model.rftn(wtype="P", nt=8, dt=0.1, output_path="bad", inca=90.0),
    "inca should be below 90 degrees",
)
expect_value_error(
    lambda: model.rftn(wtype="P", nt=8, dt=0.1, output_path="bad", rayp=0.0),
    "rayp should be positive",
)

for name in ["PY_P", "PY_S", "bad"]:
    path = Path(name)
    if path.is_dir():
        shutil.rmtree(path, ignore_errors=True)
