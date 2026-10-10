import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from obspy import read
from scipy.io import netcdf_file


for model, title in (("layered", "Layered halfspace"), ("halfspace", "Homogeneous halfspace")):
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True, layout="constrained")
    with netcdf_file(f"static_{model}.nc", mmap=False) as dataset:
        static = {component: dataset.variables[component][:].item() for component in "ZNE"}
    if model == "halfspace":
        with netcdf_file("okada.nc", mmap=False) as dataset:
            okada = {component: dataset.variables[component][:].item() for component in "ZNE"}

    for ax, component in zip(axes, "ZNE"):
        trace = read(f"syn_{model}/{component}.sac")[0]
        time = trace.times() + trace.stats.sac.b
        ax.plot(time, trace.data, color="C0", lw=1.8, label="Dynamic (greenfn + syn)")
        ax.axhline(static[component], color="C1", ls="--", lw=1.4, label="Static (static greenfn + syn)")
        print(f"{model} {component}: dynamic end = {trace.data[-1]:.6e}, static = {static[component]:.6e} cm")
        if model == "halfspace":
            lamb = read(f"lamb/{component}.sac")[0]
            ax.plot(lamb.times() + lamb.stats.sac.b, lamb.data, color="0.25", ls=":", lw=1.2, label="Dynamic (Lamb)")
            ax.axhline(okada[component], color="C3", ls="-.", lw=1.2, label="Static (Okada)")
            print(f"  Lamb end = {lamb.data[-1]:.6e}, Okada = {okada[component]:.6e} cm")

        ax.set_title(component)
        ax.set_ylabel("Displacement (cm)")
        ax.set_xmargin(0)
        ax.grid(alpha=0.3, lw=0.5)

    axes[0].legend(fontsize=8, loc="best")
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle(title)
    fig.savefig(f"compare_{model}.svg")
    plt.close(fig)
