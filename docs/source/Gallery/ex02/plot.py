import numpy as np
import matplotlib.pyplot as plt
from obspy import *
import glob, os

diff = 0.001
depsrc = 2
deprcv = 3.2
dist = 10

def plot_one(ax, fname, prefix):
    tr0 = read(path)[0]
    if prefix == 'z':
        tr1 = read(os.path.join(f"./GRN/milrow_{depsrc}_{deprcv+diff:.3f}_{dist}/", fname))[0]
    elif prefix == 'r':
        tr1 = read(os.path.join(f"./GRN/milrow_{depsrc}_{deprcv}_{dist+diff:.3f}/", fname))[0]
    else:
        return
    
    trpar = read(os.path.join(f"./GRN/milrow_{depsrc}_{deprcv}_{dist}/", prefix+fname))[0]

    sgn = -1 if prefix == 'z' else 1
    nt = tr0.stats.npts
    dt = tr0.stats.delta
    t = np.arange(nt)*dt
    ax.plot(t, (tr1.data - tr0.data) / diff * sgn, lw=0.4, c='b')
    ax.plot(t, trpar.data, lw=0.4, c='r')
    ax.set_xlim(0, nt*dt)
    ax.set_yticks([])
    ax.text(0.95, 0.9, prefix+fname, transform=ax.transAxes, ha='right', va='top')


pathLst = glob.glob(f"./GRN/milrow_{depsrc}_{deprcv}_{dist}/???.sac")
pathLst.sort()


# ====================== z ==============================
fig, axs = plt.subplots(len(pathLst), 1, figsize=(10, 1*len(pathLst)), 
                        gridspec_kw=dict(hspace=0), sharex=True)

for i, path in enumerate(pathLst):
    ax = axs[i]
    fname = os.path.basename(path)
    plot_one(ax, fname, "z")

fig.savefig("compare_z.svg", bbox_inches='tight')


# ====================== r ==============================
fig, axs = plt.subplots(len(pathLst), 1, figsize=(10, 1*len(pathLst)), 
                        gridspec_kw=dict(hspace=0), sharex=True)

for i, path in enumerate(pathLst):
    ax = axs[i]
    fname = os.path.basename(path)
    plot_one(ax, fname, "r")

fig.savefig("compare_r.svg", bbox_inches='tight')