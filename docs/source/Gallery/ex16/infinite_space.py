import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

nt = 801     # 总点数，不要求2的幂次
dt = 0.01    # 采样时间间隔(s)
t = np.arange(0, nt)*dt
rs = np.arange(0.01, 10.01, 0.01)

Vp, Vs, Rho = 5.0, 3.0, 2.2
mu = Rho * Vs**2
lam = Rho * Vp**2 - 2*mu
Delta = (lam + mu) / (lam + 3.0*mu)

depsrc = 3.    # 震源深度 km
deprcv = 1.    # 台站深度 km

def static_theoretical(ij:str):
    D = np.hypot(depsrc - deprcv, rs)
    theta = np.arctan((depsrc - deprcv) / rs)
    gamz = np.sin(theta)
    gamr = np.cos(theta)

    delta_ij = 1.0 if ij[0]==ij[1] else 0.0
    gam_ij = np.ones_like(D)
    for s in list(ij):
        if s == 'z':
            gam_ij *= gamz
        elif s == 'r':
            gam_ij *= gamr
        elif s == 't':
            gam_ij *= 0.0

    return 1/(4.0*np.pi*mu) * Delta/(1 + Delta) / D * (delta_ij / Delta + gam_ij)

def dynamic_theoretical(ij:str):
    r = 10
    R = np.hypot(depsrc - deprcv, r)
    theta = np.arctan((depsrc - deprcv) / r)
    gamz = np.sin(theta)
    gamr = np.cos(theta)

    delta_ij = 1.0 if ij[0]==ij[1] else 0.0
    gam_ij = 1.0
    for s in list(ij):
        if s == 'z':
            gam_ij *= gamz
        elif s == 'r':
            gam_ij *= gamr
        elif s == 't':
            gam_ij *= 0.0

    def heaviside(t0):
        arr = np.zeros_like(t)
        arr[t >= t0] = 1
        return arr

    # 积分后的结果
    near = 1/(4.0*np.pi*Rho) * (3.0*gam_ij - delta_ij) / R**3 * 0.5 * heaviside(R/Vp)*(np.fmin(t, R/Vs)**2 - (R/Vp)**2)
    Pfar = 1/(4.0*np.pi*Rho*Vp**2) * gam_ij / R * heaviside(R/Vp)
    Sfar = - 1/(4.0*np.pi*Rho*Vs**2) * (gam_ij - delta_ij) / R * heaviside(R/Vs)

    return near + Pfar + Sfar

# =============================================================
chLst = ['VFZ', 'VFR', 'HFZ', 'HFR', 'HFT']

staticLst = []
staticLst.append(static_theoretical('zz') * (-1))
staticLst.append(static_theoretical('zr') * (-1))
staticLst.append(static_theoretical('zr'))
staticLst.append(static_theoretical('rr'))
staticLst.append(static_theoretical('tt'))

dynamicLst = []
dynamicLst.append(dynamic_theoretical('zz') * (-1))
dynamicLst.append(dynamic_theoretical('zr') * (-1))
dynamicLst.append(dynamic_theoretical('zr'))
dynamicLst.append(dynamic_theoretical('rr'))
dynamicLst.append(dynamic_theoretical('tt'))

fig = plt.figure(figsize=(8, 6))
fig1, fig2 = fig.subfigures(1, 2, wspace=0.1)

axs1 = fig1.subplots(5, 1, sharex=True, gridspec_kw=dict(hspace=0))
axs2 = fig2.subplots(5, 1, sharex=True, gridspec_kw=dict(hspace=0))
for i in range(5):
    prop1 = dict(ls='-', c='b', lw=2)

    ax = axs1[i]
    ax.plot(t, dynamicLst[i], **prop1)
    ax.text(0.96, 0.9, chLst[i], transform=ax.transAxes, ha='right', va='top', bbox=dict(fc='w'))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax.grid()
    ax.set_xmargin(0)
    ax.set_ymargin(0.3)

    ax = axs2[i]
    ax.plot(rs, staticLst[i], **prop1)
    ax.text(0.96, 0.9, chLst[i], transform=ax.transAxes, ha='right', va='top', bbox=dict(fc='w'))
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

fig.text(0.0, 1.0, f"theoretical solution in an infinite space", 
         ha='left', va='bottom', fontsize=12, bbox=dict(fc='w'))

fig.tight_layout()
fig.savefig(f"theoretical.svg", bbox_inches='tight')
