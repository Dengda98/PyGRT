from obspy import read 
import matplotlib.pyplot as plt
import numpy as np

st = read("GRN/*/VFZ.sac")
nt = st[0].stats.npts
dt = st[0].stats.delta

scale = 0.05
t = np.arange(0, nt)*dt
fig, ax = plt.subplots(1, 1, figsize=(6, 9))
for tr in st:
    data = tr.data.copy()
    data /= np.max(np.abs(data)) # 归一化
    dist = tr.stats.sac['dist']
    
    ax.plot(data*scale+dist, t, 'k', lw=0.2) 
    # ax.fill_betweenx(t, data*scale+dist, dist, where=(data > 1e-4), color='k')

ax.set_xlabel("Distance (km)")
ax.set_ylabel("Time (s)")
ax.set_xlim([0.5, 1.5])
ax.set_ylim([0, 15])
ax.yaxis.set_inverted(True)

fig.savefig("multi_traces.png", bbox_inches='tight', dpi=100)