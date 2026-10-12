import numpy as np
from obspy import Stream, read
import matplotlib.pyplot as plt

def plot6(st6: Stream, title: str, out: str | None = None, *, zero_channels: tuple[str, ...] = ()):
    nt = st6[0].stats.npts
    dt = st6[0].stats.delta
    t = np.arange(nt)*dt

    MAX, MIN = 0, 9e30
    for tr in st6:
        d = np.abs(tr.data)
        if MAX < np.max(d):
            MAX = np.max(d)
        if MIN > np.min(d):
            MIN = np.min(d)

    travtP = st6[0].stats.sac['t0']
    travtS = st6[0].stats.sac['t1']

    fig, axs = plt.subplots(len(st6), 1, figsize=(10, 1.2*len(st6)), gridspec_kw=dict(hspace=0.0), sharex=True)
    for i in range(len(st6)):
        tr = st6[i]
        ax = axs[i]

        ax.plot(t, tr.data, c='k', lw=0.5, label=tr.stats.channel)
        ax.legend(loc='upper left')

        # 理论为零或幅值极小的分量采用整体量级设置纵轴，避免放大数值残差
        m = np.max(np.abs(tr.data))
        ylims = ax.get_ylim()
        if tr.stats.channel in zero_channels or m / MAX < 1e-5:
            ylims = np.array([-1, 1]) * MAX
            ax.set_ylim(ylims)

        # 绘制到时
        ax.vlines(travtP, *ylims, colors='b')
        ax.text(travtP, ylims[1], "P", ha='left', va='top', color='b')
        ax.vlines(travtS, *ylims, colors='r')
        ax.text(travtS, ylims[1], "S", ha='left', va='top', color='r')

        ax.margins(x=0)
        ax.set_xlim([t[0], t[-1]])
        ax.set_ylim(np.array(ylims)*1.2)

    fig.suptitle(title)

    if out is not None:
        fig.savefig(out, bbox_inches='tight')
    return fig, axs

fig, _ = plot6(read("syn_dc_zne/strain_*.sac"), "Strain", "strain.svg")
plt.close(fig)
fig, _ = plot6(read("syn_dc_zne/rotation_*.sac"), "Rotation", "rotation.svg")
plt.close(fig)
fig, _ = plot6(read("syn_dc_zne/stress_*.sac"), "Stress", "stress.svg", zero_channels=("ZZ", "ZN", "ZE"))
plt.close(fig)

st_projected = Stream()
for name, label in zip(("tau_s", "sigma_n", "coulomb"), (r"$\Delta\tau_s$", r"$\Delta\sigma_n$", r"$\Delta CFS$")):
    tr = read(f"syn_dc_zne/{name}.sac")[0]
    tr.data *= 1e-7
    tr.stats.channel = f"{label} (MPa)"
    st_projected += tr
fig, axes = plot6(st_projected, "Stress projection and Coulomb stress change")
axes[-1].set_xlabel("Time (s)")
fig.savefig("coulomb.svg", bbox_inches="tight")
plt.close(fig)
