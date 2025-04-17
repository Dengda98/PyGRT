import numpy as np
from obspy import *
import matplotlib.pyplot as plt
import pygrt

# load halfspace model
modarr = np.loadtxt("killari")
modarr2 = np.loadtxt("killari_nosoil")

depsrc = 3.0
deprcv = 0.0

rs = np.array([3])

nt = 500
dt = 0.01 

r = rs[0]
t = np.arange(0, nt)*dt

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)
pymod2 = pygrt.PyModel1D(modarr2, depsrc, deprcv)

# compute Green's Functions
st_grn = pymod.compute_grn(
    distarr=rs, 
    nt=nt, 
    dt=dt, 
    vmin_ref=3.5,
    keps=-1,
    Length=20, 
    gf_source=['DC']
)[0]
st_grn2 = pymod2.compute_grn(
    distarr=rs, 
    nt=nt, 
    dt=dt, 
    vmin_ref=3.5,
    keps=-1,
    Length=20, 
    gf_source=['DC']
)[0]

# define time function, 1/tau sin3 (\pi t/tau)
def custome_func(tau, dt, L):
    t = np.arange(0, L, dt)
    amp = 1/tau * (np.sin(np.pi*t/tau))**3
    amp[t>tau] = 0.0
    return t, amp



pygrt.utils.stream_convolve(st_grn, custome_func(0.05, dt, 0.2)[1])
pygrt.utils.stream_convolve(st_grn2, custome_func(0.05, dt, 0.2)[1])

st = pygrt.utils.gen_syn_from_gf_DC(st_grn, 1e24, 0, 45, 90, 45)
st2 = pygrt.utils.gen_syn_from_gf_DC(st_grn2, 1e24, 0, 45, 90, 45)

norm = np.linalg.norm(st.select(component='R')[0].data, np.inf)


fig, axs = plt.subplots(3, 2, figsize=(12, 6))

for i, (ch, ylabel) in enumerate(zip(['Z', 'R', 'T'], ['Vertical', 'Radial', 'Transverse'])):
    tr = st.select(component=ch)[0]
    tr2 = st2.select(component=ch)[0]
    ax = axs[i][0]
    ax2 = axs[i][1]
    ax.plot(t, tr.data/norm, 'k-', lw=1.2)
    ax2.plot(t, tr2.data/norm, 'k-', lw=1.2)

    ax.set_xlim([0, 3])
    ax2.set_xlim([0, 3])

    ax.set_ylabel(ylabel, fontsize=13)

    ax.set_ylim([-1, 1])
    ax2.set_ylim([-1, 1])

    ax.grid()
    ax2.grid()


axs[0,0].set_title("with 5m soft surface layer", fontsize=16)
axs[0,1].set_title("without the soft surface layer", fontsize=16)

axs[-1,0].set_xlabel("Time (s)", fontsize=13)
axs[-1,1].set_xlabel("Time (s)", fontsize=13)

fig.savefig('syn2_compare.pdf')