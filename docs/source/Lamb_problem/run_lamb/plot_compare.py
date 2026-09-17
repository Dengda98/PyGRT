"""比较 Lamb 解析解与 greenfn + syn 的均匀半空间动态结果"""

import shutil
from math import isfinite
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from obspy import read

# BEGIN LAMB
import pygrt

pygrt.utils.lamb(
    modelparams=(8.0, 4.62, 3.3),
    depsrc=5.0,
    deprcv=1.0,
    dist=15.0,
    nt=1300,
    dt=0.005,
    azimuth=30.0,
    output_path="lamb_python",
    scale=1e24,
    strike=33.0,
    dip=50.0,
    rake=120.0,
    time_function="t/0.1/0.1/0.2",
    zne=True,
)
# END LAMB

pymod = pygrt.PyModel1D(grn="GRN", modelpath="halfspace")
pymod.greenfn(
    depsrc=5.0,
    deprcv=1.0,
    dists=15.0,
    nt=1300,
    dt=0.005,
)
pymod.syn(
    azimuth=30.0,
    output_path="syn_python",
    scale=1e24,
    strike=33.0,
    dip=50.0,
    rake=120.0,
    time_function="t/0.1/0.1/0.2",
    zne=True,
)


def get_sac_arrivals(trace):
    """读取 SAC 头段中已经定义的震相到时"""
    arrivals = []
    sac = trace.stats.sac
    for index in range(10):
        arrival = getattr(sac, f"t{index}", None)
        phase = str(getattr(sac, f"kt{index}", "")).strip()
        if arrival is None or not isfinite(arrival) or arrival <= -12344.0:
            continue
        if not phase or phase == "-12345":
            continue
        arrivals.append((float(arrival), phase))
    return arrivals


def annotate_arrivals(axes, arrivals):
    """在所有分量上绘制震相到时，并在顶部分量标注震相名称"""
    for index, (arrival, phase) in enumerate(arrivals):
        for axis in axes:
            axis.axvline(arrival, color="red", linewidth=0.6, alpha=0.55)
        axes[0].text(
            arrival,
            1.01,
            phase,
            transform=axes[0].get_xaxis_transform(),
            color="red",
            fontsize=8,
            ha="center",
            va="bottom",
            clip_on=False,
        )


def plot_comparison(lamb_directory: Path, syn_directory: Path, output_path: Path) -> None:
    figure, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True, layout="constrained")
    arrivals = None
    for axis, component in zip(axes, "ZNE"):
        lamb = read(str(lamb_directory / f"{component}.sac"))
        syn = read(str(syn_directory / f"{component}.sac"))
        if arrivals is None:
            arrivals = get_sac_arrivals(lamb[0])
        pygrt.utils.stream_integral(lamb)
        pygrt.utils.stream_integral(syn)
        lamb = lamb[0]
        syn = syn[0]
        axis.plot(lamb.times(), lamb.data, color="0.6", label="lamb", lw=2.5)
        axis.plot(syn.times(), syn.data, "--", color="blue", label="greenfn + syn", lw=0.8)
        axis.set_ylabel(component+" (cm)")
        axis.grid(linewidth=0.2)
        axis.set_xmargin(0)

    if arrivals:
        annotate_arrivals(axes, arrivals)
    axes[-1].set_xlabel("Time (s)")
    axes[0].legend()
    figure.savefig(output_path, bbox_inches="tight")
    plt.close(figure)

plot_comparison(Path("lamb_python"), Path("syn_python"), Path("lamb_compare.svg"))

for path in (Path("GRN"), Path("lamb_python"), Path("syn_python")):
    if path.is_dir():
        shutil.rmtree(path)
