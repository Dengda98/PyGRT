import numpy as np
import matplotlib.pyplot as plt
import pygrt 

statsdata1 = pygrt.utils.read_statsfile("GRN1_grtstats/milrow_5_0/K_0500_2.50000e-01")
statsdata2 = pygrt.utils.read_statsfile("GRN2_grtstats/milrow_5_0/K_0500_2.50000e-01")
srctype = 'DD_q'

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
xlim = [0, 0.65]

karr = statsdata1['k']
_filt = np.logical_and(karr <= xlim[1], karr >= xlim[0])
karr = karr[_filt]
Farr = statsdata1[srctype][_filt]
ax1.vlines(karr, -1, 1, lw=0.1)
ax1.plot(karr, np.real(Farr)/np.max(np.real(Farr)), 'o-', ms=0.5, mec='r', mfc='r', c='k', lw=0.5)
ax1.grid()
ax1.set_xlim(xlim)
ax1.set_title(f"Fixed Interval (npts={len(karr)})")

karr = statsdata2['k']
_filt = np.logical_and(karr <= xlim[1], karr >= xlim[0])
karr = karr[_filt]
Farr = statsdata2[srctype][_filt]
ax2.vlines(karr, -1, 1, lw=0.1)
ax2.plot(karr, np.real(Farr)/np.max(np.real(Farr)), 'o-', ms=0.5, mec='r', mfc='r', c='k', lw=0.5)
ax2.grid()
ax2.set_xlim(xlim)
ax2.set_title(f"Self-Adaptive Interval (npts={len(karr)})")

fig.tight_layout()
fig.savefig("safim.png", dpi=300, bbox_inches='tight')