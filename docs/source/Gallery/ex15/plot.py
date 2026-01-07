import numpy as np
import pygrt
import matplotlib.pyplot as plt
import sys

modname = sys.argv[1]
modarr = np.loadtxt(modname)

deprcv = 0.0
depsrc = 5.0
distarr = np.arange(0.1, 50, 1)

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)

# 零频频谱, 此时 dt,zeta 变相用于控制虚频率的大小
pygrnLst, _, _ = pymod._get_grn_spectra(distarr, nt=1, dt=100, zeta=1.0, keepAllFreq=True)

# 静态解
static_grn = pymod.compute_static_grn(xarr=[0.0], yarr=distarr)

# 绘制
coef = 1e-20 * 1e25  # 1e25 为地震矩

# 绘制零频结果
fig, ax = plt.subplots(figsize=(6,6))
dynamic_Z = np.zeros_like(distarr)
dynamic_R = np.zeros_like(distarr)
for i in range(len(distarr)):
    dynamic_Z[i] = np.real(pygrnLst[i][0][0].cmplx_grn[0]) * coef * (-1)
    dynamic_R[i] = np.real(pygrnLst[i][0][1].cmplx_grn[0]) * coef

ax.plot(distarr, dynamic_Z, 'k', label='Dynamic Z')
ax.plot(distarr, dynamic_R, 'k', label='Dynamic R')

ax.plot(distarr, static_grn['EXZ'][0] * coef, 'ro', ms=4, label='Static Z')
ax.plot(distarr, static_grn['EXR'][0] * coef, 'bo', ms=4, label='Static R')

ax.set_xlim(0, 50)
ax.set_xlabel('Distance (km)')
ax.set_ylim(ymin=0)
ax.set_ylabel('Displacement (cm)')

ax.legend(loc='upper right')
ax.grid()

fig.savefig(f'{modname}.svg', bbox_inches='tight')