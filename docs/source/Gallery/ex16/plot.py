import pygrt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import sys

yarr = np.arange(0.01, 10.01, 0.01)
rs = np.array([10]) # 震中距数组，km
nt = 801     # 总点数，不要求2的幂次
dt = 0.01    # 采样时间间隔(s)
Vp, Vs, Rho = 5.0, 3.0, 2.2
t = np.arange(0, nt)*dt

bound1 = sys.argv[1]
bound2 = sys.argv[2]
isuniform = sys.argv[3] == 'y'

# =============================================================
depsrc = 3.    # 震源深度 km
deprcv = 1.    # 台站深度 km
if isuniform:
    modarr = np.array([
        [2.0,  Vp, Vs, Rho],
    ])
else:
    modarr = np.array([
        [2.0,  Vp, Vs, Rho],
        [1.6,  Vp+0.1, Vs+0.1, Rho+0.1],
        [2.8,  Vp+0.4, Vs+0.4, Rho+0.2],
        [0,    Vp+0.6, Vs+0.6, Rho+0.3],
    ])
pymod1 = pygrt.PyModel1D(modarr, depsrc, deprcv, topbound=bound1, botbound=bound2)
st1 = pymod1.compute_grn(distarr=rs, nt=nt, dt=dt, keepAllFreq=True)[0]
pygrt.utils.stream_integral(st1)
static1 = pymod1.compute_static_grn(xarr=[0.0], yarr=yarr)

# =============================================================
# 设置上下翻转模型
if bound2 == 'halfspace':
    modarr2 = np.flipud(modarr)
    H = np.sum(modarr[:-1, 0])
    modarr2[0,0] = max(max((depsrc, deprcv)) - H, 0) + 1.0  #  顶部半空间加 1km 厚度作为参考偏移
else:
    modarr2 = np.flipud(modarr[:-1, :])

H = np.sum(modarr2[:, 0])
if bound1 != 'halfspace':
    modarr2 = np.vstack((modarr2, modarr[0,:]))

depsrc2 = H - depsrc  # 震源深度 km
deprcv2 = H - deprcv  # 台站深度 km
print(modarr2)
print(modarr2.shape, depsrc2, deprcv2)

pymod2 = pygrt.PyModel1D(modarr2, depsrc2, deprcv2, topbound=bound2, botbound=bound1)  # 整理好的模型对象
st2 = pymod2.compute_grn(distarr=rs, nt=nt, dt=dt, keepAllFreq=True)[0]
static2 = pymod2.compute_static_grn(xarr=[0.0], yarr=yarr)
pygrt.utils.stream_integral(st2)


# =============================================================
chLst = ['VFZ', 'VFR', 'HFZ', 'HFR', 'HFT']

fig = plt.figure(figsize=(8, 6))
fig1, fig2 = fig.subfigures(1, 2, wspace=0.1)

axs1 = fig1.subplots(5, 1, sharex=True, gridspec_kw=dict(hspace=0))
axs2 = fig2.subplots(5, 1, sharex=True, gridspec_kw=dict(hspace=0))
for i in range(5):
    ch = chLst[i]
    sgn = -1 if ch in ['VFR', 'HFZ'] else 1

    prop1 = dict(ls='-', c='0.7', lw=3.5)
    prop2 = dict(ls='--', c='k', lw=1.2)

    ax = axs1[i]
    data = st1.select(channel=chLst[i])[0].data
    ax.plot(t, data, **prop1)
    data = st2.select(channel=chLst[i])[0].data * sgn
    ax.plot(t, data, **prop2)
    ax.text(0.96, 0.9, chLst[i], transform=ax.transAxes, 
            ha='right', va='top', bbox=dict(fc='w'))

    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.grid()
    ax.set_xmargin(0)
    ax.set_ymargin(0.3)

    ax = axs2[i]
    ax.plot(yarr, static1[chLst[i]][0], **prop1)
    ax.plot(yarr, static2[chLst[i]][0] * sgn, **prop2)
    ax.text(0.96, 0.9, chLst[i], transform=ax.transAxes, 
            ha='right', va='top', bbox=dict(fc='w'))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.grid()
    ax.set_xmargin(0)
    ax.set_ymargin(0.3)

ax = axs1[0]
ax.set_title("Dynamic displacements")
ax = axs1[-1]
ax.set_xlabel("Time (s)")

ax = axs2[0]
ax.set_title("Static displacements")
ax = axs2[-1]
ax.set_xlabel("Distance (km)")

fig.text(0.0, 1.0, f"gray-solid:  {bound1} + {bound2}\n"
                    f"black-dash:  {bound2} + {bound1}", 
         ha='left', va='bottom', fontsize=12, bbox=dict(fc='w'))

fig.tight_layout()
# fig.savefig(f"test_flip1.svg", bbox_inches='tight')
fig.savefig(f"{bound1}_{bound2}.svg", bbox_inches='tight')


