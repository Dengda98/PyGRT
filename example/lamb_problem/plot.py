from obspy import read 
import matplotlib.pyplot as plt
import numpy as np
import pygrt

st = read("GRN/*/*")
nt = st[0].stats.npts
dt = st[0].stats.delta

# convolve with step function
pygrt.utils.stream_integral(st)

t = np.arange(0, nt)*dt
fig, ax = plt.subplots(1, 1, figsize=(4, 6))
for itr, tr in enumerate(st):
    data = tr.data.copy()
    data /= np.max(np.abs(data)) # 归一化
    
    h, = ax.plot(t, data+itr*2, lw=1.5) 
    ax.text(2, itr*2, tr.stats.channel, c=h.get_color(), ha='right')

ax.set_xlabel("Time (s)")
ax.set_xlim([2.0, 2.7])
ax.set_yticks([])

for (_, spine) in ax.spines.items():
    spine.set_visible(False)

ax.spines['bottom'].set_visible(True)

fig.savefig("stream.png", bbox_inches='tight', dpi=600)