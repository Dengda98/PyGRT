import pygrt
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from obspy import read

plt.rcParams.update({
    "font.sans-serif": "Times New Roman",
    "mathtext.fontset": "cm"
})

# 定义半无限空间模型
Vp = 8     # km/s
Vs = 4.62  # km/s
Rho = 3.3  # g/cm^3

# 泊松比
nu = 0.5 * (1 - 2*(Vs/Vp)**2) / (1 - (Vs/Vp)**2)

# 模型数组，半无限空间
modarr = np.array([
    [0.,  Vp, Vs, Rho, 9e10, 9e10],
])

# 剪切模量
mu = Vs**2 * Rho

depsrc = 0.0  # 震源深度 km，位于地表
deprcv = 5.0  # 台站深度 km

rs = np.array([10]) # 震中距数组，km

nt = 1300     # 总点数，不要求2的幂次
dt = 0.005    # 采样时间间隔(s)

idx = 0      # 震中距索引
epicentral_distance = rs[idx]

r_straight = np.hypot(epicentral_distance, deprcv)
# 将频域解的物理时间转换为无量纲时间
tbar = np.arange(0, nt)*dt * Vs/r_straight

modfile = "_halfspace_mod"
np.savetxt(modfile, modarr)
pymod = pygrt.PyModel1D(grn="GRN", modelpath=modfile)
# 计算格林函数（仅一个震中距，可用通配符读回）
pymod.greenfn(depsrc=depsrc, deprcv=deprcv, dists=rs, nt=nt, dt=dt, calc_upar=True, Length=50)
st = read("GRN/*/*.sac")

# 卷积阶跃函数
pygrt.utils.stream_integral(st)


def plot(st, prefix, u, scale, sub, ylim):
    # 频域解
    v = np.zeros_like(u)
    coef = np.pi * np.pi * mu * r_straight * scale
    pref = prefix[-1] if len(prefix) > 0 else ''
    v[:,0] = st.select(channel=f'{pref}HFR')[0].data * coef
    v[:,1] = st.select(channel=f'{pref}VFR')[0].data * coef
    v[:,2] = st.select(channel=f'{pref}HFT')[0].data * coef
    v[:,3] = st.select(channel=f'{pref}HFZ')[0].data * coef * (-1)
    v[:,4] = st.select(channel=f'{pref}VFZ')[0].data * coef * (-1)

    fig, axs = plt.subplots(5, 2, figsize=(10, 10), sharex=True)
    labels = [
        rf"$\bar{{G}}^H_{{11{sub}}}$",
        rf"$\bar{{G}}^H_{{13{sub}}}$",
        rf"$\bar{{G}}^H_{{22{sub}}}$",
        rf"$\bar{{G}}^H_{{31{sub}}}$",
        rf"$\bar{{G}}^H_{{33{sub}}}$",
    ]

    for i in range(5):
        ax = axs[i, 0]
        ax.plot(tbar, u[:,i], c='0.5', lw=1.5)
        ax.set_ylim(*ylim)
        ax.set_xlim(0, 2)
        ax.grid(lw=0.4)
        ax.text(0.05, 0.92, labels[i], transform=ax.transAxes, ha='left', va='top', fontsize=12)

        ax = axs[i, 1]
        ax.plot(tbar, u[:,i], c='0.5', lw=1.5)
        ax.plot(tbar, v[:,i], c='b', lw=0.7)
        ax.set_ylim(*ylim)
        ax.set_xlim(0, 2)
        ax.grid(lw=0.4)
        ax.text(0.05, 0.92, labels[i], transform=ax.transAxes, ha='left', va='top', fontsize=12)

    axs[0,0].set_title("From Time-Domain")
    axs[0,1].set_title("From Frequency-Domain")

    fig.savefig(f"lamb2_compare_freq_time_receiver{prefix}.svg", bbox_inches='tight')




# 时域解
u, _, ur = pygrt.utils.lamb2(nu=nu, tbar=tbar, R=epicentral_distance, deprcv=deprcv, azimuth=0.0)
u = u.reshape(-1, 9)[:, [0,2,4,6,8]]
ur = ur.reshape(-1, 3, 9)[:, :, [0,2,4,6,8]]

plot(st, '', u, 1.0, "", (-2, 2))
plot(st, '_z', ur[:, 2, :], - r_straight, ",3", (-4, 4))
plot(st, '_r', ur[:, 0, :], r_straight, ",1", (-4, 4))

# 删除中间计算结果，仅保留成图
import shutil
for name in ["GRN", "_halfspace_mod"]:
    p = Path(name)
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    elif p.is_file():
        p.unlink(missing_ok=True)
