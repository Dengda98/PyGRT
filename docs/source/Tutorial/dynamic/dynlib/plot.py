import matplotlib.pyplot as plt
from obspy import read
import pygrt

# BEGIN GEOMETRY
fig, axes = pygrt.plot.plot_sCoulomb_faults_geometry("fault_rupture.inp")
for ax in axes:
    ax.get_legend().remove()
fig.savefig("fault.svg")
plt.close(fig)
# END GEOMETRY

# BEGIN WAVEFORM
fig, axes = plt.subplots(3, 1, figsize=(10, 4), gridspec_kw=dict(hspace=0.0), sharex=True)
for ax, component in zip(axes, ("Z", "N", "E")):
    tr = read(f"syn_ff/{component}.sac")[0]
    ax.plot(tr.times() + tr.stats.sac.b, tr.data, c="k", lw=0.8, label=f"{component} (cm)")
    ax.legend(loc="upper left")
    ylims = ax.get_ylim()
    ax.set_ylim([value * 1.2 for value in ylims])
    ax.margins(x=0)
fig.suptitle("Finite-fault displacement at a surface station")
axes[-1].set_xlabel("Time (s)")
fig.savefig("syn_ff.svg", bbox_inches="tight")
plt.close(fig)
# END WAVEFORM

fig, axes = plt.subplots(3, 1, figsize=(10, 4), gridspec_kw=dict(hspace=0.0), sharex=True)
for ax, name, label in zip(axes, ("tau_s", "sigma_n", "coulomb"), (r"$\Delta\tau_s$", r"$\Delta\sigma_n$", r"$\Delta CFS$")):
    tr = read(f"syn_ff/{name}.sac")[0]
    ax.plot(tr.times() + tr.stats.sac.b, tr.data * 1e-7, c="k", lw=0.8, label=f"{label} (MPa)")
    ax.legend(loc="upper left")
    ylims = ax.get_ylim()
    ax.set_ylim([value * 1.2 for value in ylims])
    ax.margins(x=0)
fig.suptitle("Finite-fault stress projection and Coulomb stress change")
axes[-1].set_xlabel("Time (s)")
fig.savefig("coulomb.svg", bbox_inches="tight")
plt.close(fig)
