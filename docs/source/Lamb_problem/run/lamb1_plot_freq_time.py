import pygrt
import numpy as np
import matplotlib.pyplot as plt

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

depsrc = 0.0  # 震源深度 km
deprcv = 0.0  # 台站深度 km

rs = np.array([10]) # 震中距数组，km

nt = 1000     # 总点数，不要求2的幂次
dt = 0.005    # 采样时间间隔(s)

idx = 0      # 震中距索引
r = rs[idx]

t = np.arange(0, nt)*dt * Vs/r 


pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)  # 整理好的模型对象
# 计算格林函数频谱
st = pymod.compute_grn(
    distarr=rs, 
    nt=nt, 
    dt=dt, 
)[0]

# 卷积阶跃函数
pygrt.utils.stream_integral(st)


# 时域解
u = pygrt.utils.solve_lamb1(nu, t, 0).reshape(-1, 9)
u = u[:, [0,2,4,6,8]]
u.shape

# 频域解
v = np.zeros_like(u)
coef = np.pi * np.pi * mu * r 
v[:,0] = st.select(channel='HFR')[0].data * coef
v[:,1] = st.select(channel='VFR')[0].data * coef
v[:,2] = st.select(channel='HFT')[0].data * coef
v[:,3] = st.select(channel='HFZ')[0].data * coef * (-1)
v[:,4] = st.select(channel='VFZ')[0].data * coef * (-1)

fig, axs = plt.subplots(5, 2, figsize=(10, 10), sharex=True)
labels = [r"$\bar{{G}}^H_{11}$", r"$\bar{{G}}^H_{13}$", r"$\bar{{G}}^H_{22}$", r"$\bar{{G}}^H_{31}$", r"$\bar{{G}}^H_{33}$"]
for i in range(5):
    ax = axs[i, 0]
    ax.plot(t, u[:,i])
    ax.set_ylim(-2, 2)
    ax.set_xlim(0, 2)
    ax.grid()
    ax.text(0.05, 0.92, labels[i], transform=ax.transAxes, ha='left', va='top', fontsize=12)

    ax = axs[i, 1]
    ax.plot(t, v[:,i])
    ax.set_ylim(-2, 2)
    ax.set_xlim(0, 2)
    ax.grid()
    ax.text(0.05, 0.92, labels[i], transform=ax.transAxes, ha='left', va='top', fontsize=12)

axs[0,0].set_title("From Time-Domain")
axs[0,1].set_title("From Frequency-Domain")

fig.savefig("lamb1_compare_freq_time.png", dpi=100, bbox_inches='tight')