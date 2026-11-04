"""绘制指定分量、指定虚拟源或全部台站的互相关距离剖面"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from obspy import read


# 输入、输出及绘图参数
INPUT          = Path("ambient_ccf")
OUTPUT         = Path("ambient_figures/profiles")
VIRTUAL_SOURCE = "STA001"                   # 虚拟源台站名；None 表示绘制全部台站
COMPONENT      = "ZZ"                       # 互相关分量，可选 "ZZ"、"RR" 或 "TT"
FIGSIZE        = (6, 6)                     # 图幅的宽和高，英寸
TRACE_SCALE    = 1.5                        # 每条互相关归一化后的绘图半振幅，km


# 台站对和距离均从 SAC 头段读取，不依赖文件名或台站排列方式
# 按距离排序后，同距的不同台站对也会逐道保留
correlations = read(str(INPUT / COMPONENT / "*.sac"))
correlations = sorted(correlations, key=lambda trace: trace.stats.sac.dist)

names = set()
for trace in correlations:
    names.add(trace.stats.sac.kevnm.strip())
    names.add(trace.stats.station)

sources = sorted(names) if VIRTUAL_SOURCE is None else [VIRTUAL_SOURCE]
ymax = max(trace.stats.sac.dist for trace in correlations) + TRACE_SCALE + 0.5

OUTPUT.mkdir(parents=True, exist_ok=True)


for name in sources:
    fig, ax = plt.subplots(figsize=FIGSIZE)

    for trace in correlations:
        first  = trace.stats.sac.kevnm.strip()
        second = trace.stats.station

        # 只取包含当前虚拟源的台站对
        if name not in (first, second):
            continue

        lags = trace.times() + trace.stats.sac.b
        values = trace.data.astype(float)

        # 保存方向与当前虚拟源方向相反时，仅在绘图时翻转数据
        if name == second:
            values = values[::-1]

        # 每道独立归一化，再按实际台间距平移，纵轴单位为 km
        values /= np.max(np.abs(values))
        baseline = trace.stats.sac.dist
        wiggle = baseline + TRACE_SCALE * values

        ax.plot(lags, wiggle, color="black", linewidth=0.6)
        ax.fill_between(lags, baseline, wiggle, where=values > 0, color="black", alpha=0.35, linewidth=0)

    ax.set(title=COMPONENT, xlabel="Lag (s)", ylabel="Epicentral distance (km)", xlim=(lags[0], lags[-1]), ylim=(0, ymax))

    fig.suptitle(f"Virtual source: {name}")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUTPUT / f"{name}.png", dpi=180)
    plt.close(fig)

    print(f"Wrote {OUTPUT / f'{name}.png'}")
