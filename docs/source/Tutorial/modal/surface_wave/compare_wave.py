import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from obspy import *
import sys

modetype = sys.argv[1]  # 0, 1, 2, all

st = read("GRN/*/[DS][DS]*")
st_nm = read(f"GRN_NM_{modetype}/*/[DS][DS]*")

tr = st[0]
nt = tr.stats.npts
dt = tr.stats.delta
t = np.arange(0, nt)*dt

fig = plt.figure(figsize=(9, 4), layout='constrained')
gs = GridSpec(4, 3, hspace=0, wspace=0.1, height_ratios=[1,1,1,0.3], figure=fig)
for isrc, (stype, sname) in enumerate(zip(['DD', 'DS', 'SS'], ['45°DS', 'DS', 'SS'])):
    for ich, (ch, chname) in enumerate(zip(['Z', 'R', 'T'], ['Vertical', 'Radial', 'Transverse'])):
        ax = fig.add_subplot(gs[ich, isrc])
        # 去除所有边框
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines:
            ax.spines[spine].set_visible(False)

        # 标注分量名
        if isrc == 0:
            ax.text(-0.1, 0.5, chname, va='center', ha='right', fontsize=10, transform=ax.transAxes)

        # 标注震源名
        if ich == 0:
            ax.text(0.5, 1.2, sname, va='bottom', ha='center', fontsize=10, transform=ax.transAxes)

        try:
            tr = st.select(channel=f"{stype}{ch}")[0]
        except:
            continue
        tr_nm = st_nm.select(channel=tr.stats.channel)[0]
        ax.plot(t, tr.data, 'k--', lw=0.3, label=f"DWM_{tr.stats.channel}")
        ax.plot(t, tr_nm.data, 'k-', lw=0.6, label=f"NMM_{tr.stats.channel}")

        ax.set_xlim([20, 80])
        

ax = fig.add_subplot(gs[-1, -2])
# 去除所有边框
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines:
    ax.spines[spine].set_visible(False)
ax.set_xlim([20, 80])
ax.spines['bottom'].set_visible(True)
ax.set_xticks(np.arange(20, 90, 20))
ax.set_xlabel('Time (s)')

fig.savefig(f"compare_{modetype}.svg")