"""对台站 ZNE 噪声去均值、去趋势和带通滤波，再按固定长度切段"""

import shutil
from pathlib import Path

from obspy import Stream, Trace, read


# 输入、临时目录及预处理和分段参数
INPUT          = Path("ambient_run")
OUTPUT         = Path("ambient_segments")
FREQMIN        = 0.1                        # 带通滤波频带下限，Hz
FREQMAX        = 1.0                        # 带通滤波频带上限，Hz
WINDOW_SECONDS = 300.0                      # 不重叠的分段长度，s


station_files = sorted((INPUT / "noise").glob("*/Z.sac"))


# 每次切段重建输出目录，分段文件保留供后续重复处理
if OUTPUT.exists():
    shutil.rmtree(OUTPUT)

OUTPUT.mkdir()


# 逐台站读取三分量完整记录，切段过程不依赖台站的位置和排列方式
for station_file in station_files:
    station_dir = station_file.parent
    name = station_dir.name

    stream = Stream()
    for component in "ZNE":
        stream += read(str(station_dir / f"{component}.sac"))

    # 如需降采样，可在去均值和去趋势之前进行
    # stream.decimate(...)

    # 按处理后的采样间隔和记录长度确定分段，末尾不足一窗的部分舍弃
    dt      = stream[0].stats.delta
    npts    = round(WINDOW_SECONDS / dt)
    nwindow = stream[0].stats.npts // npts

    # 在完整记录上去均值和线性趋势
    stream.detrend("demean")
    stream.detrend("linear")

    # 在完整记录上做一次零相位带通，再切段，避免每个短窗分别滤波产生边缘效应
    stream.filter("bandpass", freqmin=FREQMIN, freqmax=FREQMAX, corners=4, zerophase=True)

    for index in range(nwindow):
        first  = index * npts
        last   = first + npts
        window = OUTPUT / f"{index + 1:06d}"
        window.mkdir(exist_ok=True)

        for trace in stream:
            segment = Trace(header=trace.stats.copy())
            segment.data = trace.data[first:last].copy()
            segment.stats.starttime = trace.stats.starttime + first * dt

            segment.write(str(window / f"{name}.{trace.stats.channel}.sac"), format="SAC")

    print(f"{name}: wrote {nwindow} ZNE windows to {OUTPUT}")
