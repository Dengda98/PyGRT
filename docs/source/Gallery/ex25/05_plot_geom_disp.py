"""仅使用 synnoise 输出，绘制源台空间分布和相速度频散"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import pygrt


# 输入、输出及图幅参数
INPUT              = Path("ambient_run")
OUTPUT             = Path("ambient_figures")
GEOMETRY_FIGSIZE   = (9, 8)                 # 源台分布图的宽和高，英寸
DISPERSION_FIGSIZE = (12, 5)                # 相速度频散图的宽和高，英寸
STATION_SIZE       = 30                     # 台站标记大小


OUTPUT.mkdir(exist_ok=True)

# 从输出文件读取实际源台坐标，不预设源区形状或台阵走向
sources = np.loadtxt(INPUT / "sources.info", skiprows=1, usecols=(1, 2), ndmin=2)

# 台站文件包含采样参数头和坐标，只将坐标部分交给 loadtxt
# 空行和注释由 loadtxt 自行跳过
rows = []
with (INPUT / "receivers.info").open() as stream:
    for line in stream:
        fields = line.split()
        if fields and fields[0] in ("DT", "LENGTH"):
            continue

        rows.append(line)

stations = np.loadtxt(rows, usecols=(1, 2), ndmin=2)


# 源台空间分布，台站之间不绘制连线，也不标注台站名
fig, ax = plt.subplots(figsize=GEOMETRY_FIGSIZE)

ax.scatter(sources[:, 1], sources[:, 0], s=2, alpha=0.4, color="#165a96", label="Sources")
ax.scatter(stations[:, 1], stations[:, 0], marker="^", s=STATION_SIZE, color="black", label="Stations")

ax.set(xlabel="East (km)", ylabel="North (km)", title="Source and station distribution")
ax.set_aspect("equal")
ax.legend()

fig.tight_layout()
fig.savefig(OUTPUT / "geometry.png", dpi=180)
plt.close(fig)


# 直接使用 synnoise 输出的相速度频散文件
fig, axes = plt.subplots(1, 2, figsize=DISPERSION_FIGSIZE, sharey=True)

for ax, wave, title in zip(axes, "RL", ("Rayleigh waves", "Love waves")):
    dispersion = pygrt.io.read_dispersion(INPUT / f"phase{wave}.nc")

    # 按文件中实际包含的阶数绘制所有模态
    for mode, velocity in dispersion["c"].items():
        ax.plot(dispersion["freq"][mode], velocity, color="black", linewidth=1.2)

    ax.set(xlabel="Frequency (Hz)", title=title)
    ax.grid(alpha=0.2)
    ax.set_xmargin(0)

axes[0].set_ylabel("Phase velocity (km/s)")

fig.tight_layout()
fig.savefig(OUTPUT / "dispersion.png", dpi=180)
plt.close(fig)
