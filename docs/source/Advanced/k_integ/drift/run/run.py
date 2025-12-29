# BEGIN 1
import numpy as np
import pygrt

modarr = np.loadtxt("milrow")

pymod = pygrt.PyModel1D(modarr, depsrc=10, deprcv=0.0)

nt = 500
dt = 10

st_grn = pymod.compute_grn(5000, nt=nt, dt=dt, keepAllFreq=True, statsfile="pygrtstats")[0]
# END 1

# 仅绘制一个分量做示例
CHANNEL = "SSR"

# =================================================================
# 绘制波形
import matplotlib.pyplot as plt
from obspy import *
plt.rcParams.update({
    "font.sans-serif": "Times New Roman",
    "mathtext.fontset": "cm"
})

t = np.arange(nt)*dt
freqs = np.arange(nt//2 + 1)/(nt*dt)

def plot_waves(tr:Trace):
    fig, ax = plt.subplots(figsize=(8, 1.5))
    ax.plot(t, tr.data, 'k-', lw=0.5)
    ax.grid()
    ax.text(0.02, 0.9, CHANNEL, transform=ax.transAxes, fontsize=10, ha='left', va='top', bbox=dict(lw=1, fc='w'))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.set_xlim(0, nt*dt)
    ax.set_xlabel("Time (s)", fontsize=12)
    return fig, ax

fig, ax = plot_waves(st_grn.select(channel=CHANNEL)[0])
fig.savefig("grn.svg", bbox_inches='tight')
# =================================================================

# =================================================================
# 回退虚频率补偿
st_grn2 = st_grn.copy()
wI = st_grn2[0].stats.sac['user0']
for tr in st_grn2:
    tr.data[:] *= np.exp(-t*wI)

fig, ax = plot_waves(st_grn2.select(channel=CHANNEL)[0])
fig.savefig("grn2.svg", bbox_inches='tight')
# =================================================================


# =================================================================
# 绘制频谱
from scipy.fft import rfft, irfft

def plot_freqs(tr:Trace):
    fig, ax = plt.subplots(figsize=(8, 1.5))
    data = tr.data.copy()

    freqs0 = freqs.copy()
    freqs0[0] *= 0.1
    ax.plot(freqs0[:-1], np.abs(rfft(data))[:-1], 'k-', lw=0.8)
    ax.grid()
    ax.text(0.02, 0.9, CHANNEL, transform=ax.transAxes, fontsize=10, ha='left', va='top', bbox=dict(lw=1, fc='w'))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-4, 0.5/dt)
    ax.set_xlabel("Frequency (Hz)", fontsize=12)

    return fig, ax

fig, ax = plot_freqs(st_grn.select(channel=CHANNEL)[0])
fig.savefig("grn_freqs.svg", bbox_inches='tight')
fig, ax = plot_freqs(st_grn2.select(channel=CHANNEL)[0])
fig.savefig("grn2_freqs.svg", bbox_inches='tight')
# =================================================================


# =================================================================
# 读入核函数
import glob
paths = glob.glob("pygrtstats/K_000[0-5]_*")
paths.sort()
print(paths)

fig, axs = plt.subplots(len(paths), 1, figsize=(8, len(paths)*1), gridspec_kw=dict(hspace=0), sharex=True)
for i in range(len(paths)):
    ax = axs[i]
    path = paths[i]
    freq = float(path.split("_")[-1])
    D = pygrt.utils.read_statsfile(path)
    ax.plot(D['k'], np.real(D['SS_q']), 'k-', lw=1, label='Real')
    ax.plot(D['k'], np.imag(D['SS_q']), 'k--', lw=1, label='Imag')
    ax.text(0.98, 0.9, f"{freq:.1e} Hz", transform=ax.transAxes, fontsize=10, ha='right', va='top', bbox=dict(lw=1, fc='w'))
    ax.grid()
    ax.set_xlim(xmax=D['k'][-1])
    if i == 0:
        ax.legend(loc='lower center', ncols=2, bbox_to_anchor=(0.5, 1.05))

fig.savefig("kernels.svg", bbox_inches='tight')
# =================================================================



# =================================================================
# 跳过频段，重新计算
st_grn3 = pymod.compute_grn(5000, nt=nt, dt=dt)[0]

srctypes = ['EX', 'VF', 'HF', 'DD', 'DS', 'SS']
chlst = ['Z', 'R', 'T']
def plot_all_waves(st_grn:Stream):

    fig = plt.figure(figsize=(14, 2*len(srctypes)))
    gs = fig.add_gridspec(len(srctypes), 3, hspace=0.3, wspace=0.1)
    for ik, tr in enumerate(st_grn):
        channel = tr.stats.channel
        srctype = channel[:2]
        ch = channel[-1]
        
        ax = fig.add_subplot(gs[srctypes.index(srctype), chlst.index(ch)])

        data = tr.data.copy()

        ax.plot(t, data, 'k-', lw=0.6)
        ax.grid()
        ax.text(0.1, 0.9, channel, transform=ax.transAxes, ha='center', va='top', bbox=dict(lw=1, fc='w'))
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    
    return fig

fig = plot_all_waves(st_grn3.copy())
fig.savefig("grn3.svg", bbox_inches='tight')
# =================================================================
