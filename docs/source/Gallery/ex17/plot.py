import numpy as np
from obspy import *
import matplotlib.pyplot as plt

import pygrt

st = read("GRN/*/HF?.sac")
st_Q = read("GRN_Q/*/HF?.sac")

nt = st[0].stats.npts
dt = st[0].stats.delta

for _ in range(1):
    st.integrate()
    st_Q.integrate()



ts0 = 5 / 4.62
t = np.arange(nt)*dt / ts0

fig, axs = plt.subplots(2, 3, figsize=(12, 6), sharex=True, gridspec_kw=dict(hspace=0))

def _plot(ax, st, channel):
    tr = st.select(channel=channel)[0]
    tp = tr.stats.sac['t0'] / ts0
    ts = 1.0 
    ax.plot(t, tr.data, c='k', lw=1.0)
    ylim = ax.get_ylim()
    ax.vlines(tp, *ylim, colors='b', ls='--', lw=0.5)
    ax.vlines(ts, *ylim, colors='r', ls='--', lw=0.5)
    ax.text(0.05, 0.9, channel, transform=ax.transAxes, bbox=dict(fc='w'), ha='left')
    ax.set_xmargin(0)
    ax.ticklabel_format(axis='y', style='sci', scilimits=[0,0])

_plot(axs[0,0], st,   'HFZ')
_plot(axs[1,0], st_Q, 'HFZ')
_plot(axs[0,1], st,   'HFR')
_plot(axs[1,1], st_Q, 'HFR')
_plot(axs[0,2], st,   'HFT')
_plot(axs[1,2], st_Q, 'HFT')

axs[0,0].set_ylabel("No attenuation")
axs[1,0].set_ylabel(r"$Q_P=80, Q_S=20$")

fig.savefig("compare.svg", bbox_inches='tight')