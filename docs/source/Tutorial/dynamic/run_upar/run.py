# BEGIN GRN
import numpy as np
import pygrt 

modarr = np.loadtxt("milrow")

pymod = pygrt.PyModel1D(modarr, depsrc=2.0, deprcv=0.0)

# 传入calc_upar=True计算空间导数
stgrn = pymod.compute_grn(distarr=[10], nt=500, dt=0.02, calc_upar=True)[0]
print(stgrn.__str__(extended=True))
# 45 Trace(s) in Stream:
# .SYN..EXZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..VFZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DDZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..HFZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DSZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..SSZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zEXZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rEXZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zVFZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rVFZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDDZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rDDZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zHFZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rHFZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDSZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rDSZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zSSZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rSSZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..EXR  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..VFR  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DDR  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..HFR  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DSR  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..SSR  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zEXR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rEXR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zVFR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rVFR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDDR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rDDR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zHFR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rHFR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDSR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rDSR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zSSR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rSSR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..HFT  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DST  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..SST  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zHFT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rHFT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDST | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rDST | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zSST | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rSST | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# END GRN



# BEGIN SYN DC
# 传入calc_upar=True计算空间导数
stsyn = pygrt.utils.gen_syn_from_gf_DC(stgrn, M0=1e24, strike=33, dip=50, rake=120, az=30, calc_upar=True)
print(stsyn)
# 12 Trace(s) in Stream:
# .SYN..DCZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DCR  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DCT  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDCZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDCR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDCT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rDCZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rDCR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..rDCT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..tDCZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..tDCR | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..tDCT | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# END SYN DC


# BEGIN ZNE
# 传入ZNE=True可返回ZNE分量
stsyn = pygrt.utils.gen_syn_from_gf_DC(stgrn, M0=1e24, strike=33, dip=50, rake=120, az=30, calc_upar=True, ZNE=True)
print(stsyn)
# 12 Trace(s) in Stream:
# .SYN..DCZ  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DCN  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..DCE  | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDCZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDCN | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..zDCE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..nDCZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..nDCN | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..nDCE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..eDCZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..eDCN | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..eDCE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# END ZNE

from obspy import *
import matplotlib.pyplot as plt

def plot6(st6:Stream, title:str, out:str|None=None):
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

    travtP = stsyn[0].stats.sac['t0']
    travtS = stsyn[0].stats.sac['t1']

    fig, axs = plt.subplots(6, 1, figsize=(10, 7), gridspec_kw=dict(hspace=0.0), sharex=True)
    for i in range(6):
        tr = st6[i]
        ax = axs[i]

        ax.plot(t, tr.data, c='k', lw=0.5, label=tr.stats.channel)
        ax.legend(loc='upper left')

        ylims = ax.get_ylim()
        if np.max(np.abs(np.array(ylims)))/MAX < 1e-5:
            ylims = [-1, 1]
            ax.set_ylim(ylims)
            
        # 绘制到时
        ax.vlines(travtP, *ylims, colors='b')
        ax.text(travtP, ylims[1], "P", ha='left', va='top', color='b')
        ax.vlines(travtS, *ylims, colors='r')
        ax.text(travtS, ylims[1], "S", ha='left', va='top', color='r')

        ax.set_xlim([t[0], t[-1]])
        ax.set_ylim(np.array(ylims)*1.2)

    axs[0].set_title(title)

    if out is not None:
        fig.savefig(out, dpi=100)
    


# BEGIN STRAIN
st_strain = pygrt.utils.compute_strain(stsyn)
print(st_strain)
# 6 Trace(s) in Stream:
# .SYN..ZZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..ZN | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..ZE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..NN | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..NE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..EE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# END STRAIN
plot6(st_strain, "Strain", "strain.png")

# BEGIN STRESS
st_stress = pygrt.utils.compute_stress(stsyn)
print(st_stress)
# 6 Trace(s) in Stream:
# .SYN..ZZ | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..ZN | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..ZE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..NN | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..NE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# .SYN..EE | 1970-01-01T00:00:00.000000Z - 1970-01-01T00:00:09.980000Z | 50.0 Hz, 500 samples
# END STRESS
plot6(st_stress, "Stress", "stress.png")
