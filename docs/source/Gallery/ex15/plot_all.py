import numpy as np
import pygrt
import matplotlib.pyplot as plt
import sys

modname = sys.argv[1]
modarr = np.loadtxt(modname)

deprcv = 0.0
depsrc = 5.0
distarr = np.arange(0.1, 50, 0.5)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)

# 零频频谱, 此时 dt,zeta 变相用于控制虚频率的大小
pygrnLst, _, _ = pymod._get_grn_spectra(distarr, nt=1, dt=500, zeta=1.0, keepAllFreq=True)

# 静态解
static_grn = pymod.compute_static_grn(xarr=[0.0], yarr=distarr)

# 绘制零频结果
fig, axs = plt.subplots(2, 3, figsize=(12, 8), gridspec_kw=dict(hspace=0.3, wspace=0.3)) # 
axs = axs.ravel()

for isrc, (src, src2) in enumerate(zip(['EX', 'VF', 'HF', 'DD', 'DS', 'SS'],
                    ['Explosion', 'Vertical Force', 'Horizontal Force', '45°-Dip Slip', '90°-Dip Slip', 'Strike Slip'])):
    ax = axs[isrc]

    if src in ['VF', 'HF']:
        coef = 1e-15 * 1e20  # 1e20 dyne
    else:
        coef = 1e-20 * 1e25  # 1e25 dyne·cm

    dynamic_Z = np.zeros_like(distarr)
    dynamic_R = np.zeros_like(distarr)
    dynamic_T = np.zeros_like(distarr)
    for i in range(len(distarr)):
        dynamic_Z[i] = np.real(pygrnLst[i][isrc][0].cmplx_grn[0]) * coef * (-1)
        dynamic_R[i] = np.real(pygrnLst[i][isrc][1].cmplx_grn[0]) * coef
        dynamic_T[i] = np.real(pygrnLst[i][isrc][2].cmplx_grn[0]) * coef

    ms = 2
    lw = 0.8
    ax.plot(distarr, dynamic_Z, 'k', lw=lw, label='Dynamic Z')
    ax.plot(distarr, static_grn[f'{src}Z'][0] * coef, 'ro', ms=ms, label='Static Z')
    ax.plot(distarr, dynamic_R, 'k', lw=lw, label='Dynamic R')
    ax.plot(distarr, static_grn[f'{src}R'][0] * coef, 'bo', ms=ms, label='Static R')
    if src not in ['EX', 'VF', 'DD']:
        ax.plot(distarr, dynamic_T, 'k', lw=lw, label='Dynamic T')
        ax.plot(distarr, static_grn[f'{src}T'][0] * coef, 'go', ms=ms, label='Static T')

    ax.set_xlim(0, 50)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Displacement (cm)')
    ax.set_title(src2)
    ax.grid()

ax = axs[2]
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 0.92), 
          ncols=3, fontsize='large', columnspacing=4.0,
          bbox_transform=fig.transFigure)


fig.savefig(f'{modname}_all.svg', bbox_inches='tight')