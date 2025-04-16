import numpy as np
import matplotlib.pyplot as plt
from obspy import *
import pygrt


# load model
modarr = np.loadtxt("milrow")

depsrc = 2.0   
deprcv = 0.0  
pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)   

rs = np.array([1800]) 

nt = 1400     
dt = 1        

L1 = 20       
L2 = -5  
T0 = 190     

st_grn1 = pymod.compute_grn(
    distarr=rs, 
    nt=nt, 
    dt=dt, 
    freqband=[5/(nt*dt),-1],
    Length=L1, 
    delayT0=T0
)[0]

st_grn2 = pymod.compute_grn(
    distarr=rs, 
    nt=nt, 
    dt=dt, 
    freqband=[5/(nt*dt),-1],
    Length=L2, 
    delayT0=T0
)[0]

r = rs[0]

# integral 
_ = pygrt.utils.stream_integral(st_grn1)
_ = pygrt.utils.stream_integral(st_grn2)


t = np.arange(0, nt)*dt + T0
fig, ax = plt.subplots(1, 1, figsize=(12, 5), gridspec_kw=dict(hspace=0),sharex=True)

offset = 1.8
xlim = [200, 1400]
data1 = st_grn1.select(channel='EXZ')[0].data
norm = np.linalg.norm(data1, np.inf)
ax.plot(t, data1/norm + offset, 'b', lw=0.5, label=f'NO FIM')
ax.text(xlim[0] + (xlim[1]-xlim[0])*0.03, data1[-1]/norm + offset + 0.2, f"DWM(L={L1}r)", c='b', ha='left', va='bottom', fontsize=15)

data2 = st_grn2.select(channel='EXZ')[0].data
ax.plot(t, data2/norm, 'r', lw=0.5, label=f'FIM')
ax.text(xlim[0] + (xlim[1]-xlim[0])*0.03, data2[-1]/norm + 0.2, f"FIM(L={-L2}r)", c='r', ha='left', va='bottom', fontsize=15)

ax.plot(t, (data1-data2)/norm - offset, 'k', lw=0.5)
ax.text(xlim[0] + (xlim[1]-xlim[0])*0.03, (data1-data2)[0]/norm - offset + 0.2, f"Absolute Error", c='k', ha='left', va='bottom', fontsize=15)

ax.set_xlim(xlim)
ax.set_yticks([])
# ax.legend(fontsize=15)
ax.tick_params(labelsize=13)
ax.set_xlabel('Time (s)', fontsize=13)

ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.savefig(f"FIM_compare2.pdf")