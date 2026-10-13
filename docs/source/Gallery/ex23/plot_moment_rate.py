"""绘制 PyGRT 总标量矩率与 USGS 发布曲线的对比图"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from obspy import read


# 两条曲线分别使用各自的时间采样，USGS 文件首行给出的间隔为 0.01 s
reference_time, reference_rate = np.loadtxt("moment_rate.mr", skiprows=2).T
signal = read("syn/sig.sac")[0]
syn_time = signal.times() + signal.stats.sac.b
syn_rate = signal.data.astype(float) * 1e-7  # dyne cm/s 转为 N m/s

# 采用 4 s 的等效矩形窗展示子断层内部传播平滑，该窗宽不是 PARAM 文件给出的参数
# 离散卷积权重之和为 1，使用完整卷积保留尾部，保证总地震矩不变
box_duration = 4.0
box_samples = max(1, int(round(box_duration / signal.stats.delta)))
boxcar = np.ones(box_samples) / box_samples
smoothed_rate = np.convolve(syn_rate, boxcar, mode="full")
smoothed_time = signal.stats.sac.b + np.arange(len(smoothed_rate)) * signal.stats.delta

plt.rcParams.update({"font.size": 11, "svg.fonttype": "none"})
fig, ax = plt.subplots(figsize=(9, 4.6), layout="constrained")
ax.plot(reference_time, reference_rate / 1e19, color="0.25", ls="--", lw=3.0, label="USGS published")
ax.plot(syn_time, syn_rate / 1e19, color="C1", lw=1.6, label="PyGRT syn")
ax.plot(smoothed_time, smoothed_rate / 1e19, color="C2", lw=1.6, label="PyGRT syn + boxcar")
ax.set_xlim(0, 60)
ax.set_ylim(bottom=0)
ax.set_xlabel("Time since rupture onset (s)")
ax.set_ylabel(r"Moment rate ($10^{19}$ Nm/s)")
ax.grid(alpha=0.25, lw=0.5)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="upper left", frameon=False)
fig.savefig("moment_rate.svg")
plt.close(fig)
