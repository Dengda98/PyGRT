import matplotlib.pyplot as plt
import numpy as np
import pygrt
from obspy import read

# -----------------------------------------------------------------------------------
# BEGIN PYTHON RCVFN
model = pygrt.PyModel1D(modelpath="mod1")

model.rcvfn(
    wtype="P",
    rayp=0.03,
    nt=500,
    dt=0.1,
    output_path="PY_P",
    alp=5.0,
    delay=-10.0,
    write_components=True,
    upsampling_n=5,
)

model.rcvfn(
    wtype="S",
    inca=10.0,
    idx=0,
    nt=500,
    dt=0.1,
    output_path="PY_S",
    alp=5.0,
    delay=-10.0,
    write_components=True,
    upsampling_n=5,
)
# END PYTHON RCVFN
# -----------------------------------------------------------------------------------


# Plotting-only settings
ARRIVAL_TIME = 0.0
PLOT_SPAN = 12.0
OFFSET = 5.0
RCVFN_COLOR = "tab:red"
COMPONENT_COLOR = "tab:blue"
PHASE_LINE_COLOR = "0.45"
PLOT_CASES = (
    {
        "prefix": "P",
        "title": "P incidence",
        "xlim": (ARRIVAL_TIME - OFFSET, ARRIVAL_TIME + PLOT_SPAN),
        "xlabel": "Time after P arrival (s)",
        "phase_label_heights": (0.97, 0.82),
        "phase_label_va": "top",
    },
    {
        "prefix": "S",
        "title": "SV incidence",
        "xlim": (ARRIVAL_TIME - PLOT_SPAN, ARRIVAL_TIME + OFFSET),
        "xlabel": "Time before S arrival (s)",
        "phase_label_heights": (0.20, 0.05),
        "phase_label_va": "bottom",
    },
)


def read_model(path: str):
    model = np.loadtxt(path, comments="#", usecols=(0, 1, 2), ndmin=2)
    if model.shape[0] < 2 or model[-1, 0] != 0.0:
        raise ValueError("The model must end with a zero-thickness halfspace")
    thickness = model[:-1, 0]
    vp = model[:-1, 1]
    vs = model[:-1, 2]
    depths = np.cumsum(thickness)
    return thickness, vp, vs, depths


def vertical_travel_times(velocity: np.ndarray, rayp: float, thickness: np.ndarray):
    slowness_square = 1.0 / velocity**2 - rayp**2
    if np.any(slowness_square < 0.0):
        raise ValueError(f"Ray parameter {rayp:g} is not propagating in the model")
    return thickness * np.sqrt(slowness_square)


def phase_arrivals(case: dict):
    prefix = case["prefix"]
    trace = read(f"PY_{prefix}/{prefix}_rcvfn.sac")[0]
    rayp = float(trace.stats.sac.user1)
    thickness, vp, vs, depths = read_model("mod1")
    p_times = vertical_travel_times(vp, rayp, thickness)
    s_times = vertical_travel_times(vs, rayp, thickness)
    direct_name = "P" if prefix == "P" else "S"
    converted_name = "Ps" if prefix == "P" else "Sp"
    direct_times = p_times if prefix == "P" else s_times
    phases = [{"name": direct_name, "time": ARRIVAL_TIME, "depth": None}]

    for i, depth in enumerate(depths):
        above = slice(0, i + 1)
        below = slice(i + 1, None)
        if prefix == "P":
            time = p_times[below].sum() + s_times[above].sum() - direct_times.sum()
        else:
            time = s_times[below].sum() + p_times[above].sum() - direct_times.sum()
        phases.append({"name": converted_name, "time": time, "depth": depth})
    return phases


def read_trace(directory: str, name: str):
    trace = read(f"{directory}/{name}")[0]
    time = trace.stats.sac.b + np.arange(trace.stats.npts) * trace.stats.delta
    data = trace.data.astype(float)
    scale = np.max(np.abs(data))
    if scale == 0.0 or not np.all(np.isfinite(data)):
        raise ValueError(f"Invalid SAC data in {directory}/{name}")
    return time, data / scale


rows = (
    ("{}_rcvfn.sac", "Receiver function", RCVFN_COLOR),
    ("{}_Z.sac", "Vertical response", COMPONENT_COLOR),
    ("{}_R.sac", "Radial response", COMPONENT_COLOR),
)


def annotate_phases(ax, case: dict) -> None:
    xmin, xmax = case["xlim"]
    for iphase, phase in enumerate(case["phases"]):
        time = phase["time"]
        if not xmin <= time <= xmax:
            continue
        label = phase["name"]
        if phase["depth"] is not None:
            label += f"\n{phase['depth']:g} km"
        ax.axvline(time, color=PHASE_LINE_COLOR, linewidth=0.7, linestyle="--")
        label_height = case["phase_label_heights"][iphase % 2]
        ax.text(
            time, label_height, label,
            transform=ax.get_xaxis_transform(),
            ha="center", va=case["phase_label_va"], fontsize=8, color=PHASE_LINE_COLOR,
            bbox={"facecolor": "white", "alpha": 0.65, "edgecolor": "none", "pad": 0.2},
        )


def plot_case(case: dict) -> None:
    prefix = case["prefix"]
    fig, axes = plt.subplots(
        len(rows), 1,
        figsize=(8.0, 5.0),
        sharex=True,
        gridspec_kw={"hspace": 0.0},
    )
    for irow, (pattern, ylabel, color) in enumerate(rows):
        time, data = read_trace(f"PY_{prefix}", pattern.format(prefix))
        ax = axes[irow]
        ax.plot(time, data, color=color, linewidth=0.8)
        ax.axhline(0.0, color="0.7", linewidth=0.6)
        if irow == 0:
            annotate_phases(ax, case)
        else:
            ax.axvline(ARRIVAL_TIME, color="0.7", linewidth=0.6, linestyle=":")
        ax.set_xlim(*case["xlim"])
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color="0.9", linewidth=0.5)
        if irow == 0:
            ax.set_title(case["title"], fontsize=10)
        if irow == len(rows) - 1:
            ax.set_xlabel(case["xlabel"])

    fig.savefig(f"{prefix}_rcvfn.svg", bbox_inches="tight")
    plt.close(fig)


for case in PLOT_CASES:
    case["phases"] = phase_arrivals(case)


for case in PLOT_CASES:
    plot_case(case)
