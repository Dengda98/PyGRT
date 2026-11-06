"""读取 ZNE 分段，逐段 one-bit 和白化，在频域计算并叠加互相关"""

import shutil
from pathlib import Path

import numpy as np
from obspy import Trace, read
from scipy.fft import next_fast_len
from scipy.signal.windows import tukey


# 输入、输出及互相关参数
INPUT          = Path("ambient_segments")
OUTPUT         = Path("ambient_ccf")
FREQMIN        = 0.05                       # 白化频带下限，Hz
FREQMAX        = 2.0                        # 白化频带上限，Hz
MAX_LAG        = 50.0                       # 保留的最大正负时延，s
TAPER_FRACTION = 0.05                       # 白化频带两端余弦渐变各占的比例


def read_segment(directory, name):
    record = {}
    for component in "ZNE":
        trace = read(str(directory / f"{name}.{component}.sac"))[0]
        record[component] = trace.data.astype(float)

    return record


def process(data, nfft, taper):
    values    = np.sign(data)  # one-bit
    
    spectrum  = np.fft.rfft(values, n=nfft)
    amplitude = np.abs(spectrum)

    # 对非零幅值频点完全白化，零幅值频点保持为零
    np.divide(spectrum, amplitude, out=spectrum, where=amplitude > 0)

    # 频带外置零，白化频带两端做余弦渐变，中间频点保持完全白化
    spectrum *= taper

    return spectrum


def rotate_horizontal(north, east, azimuth):
    cosine, sine = np.cos(azimuth), np.sin(azimuth)

    # 切向为径向顺时针旋转 90 度
    radial     =  north * cosine + east * sine
    transverse = -north * sine   + east * cosine

    return radial, transverse


# 台站位置、标识和采样信息均从首个时间窗的 SAC 头段读取
windows = sorted(path for path in INPUT.glob("*") if path.is_dir())
headers = [read(str(path), headonly=True)[0].stats for path in sorted(windows[0].glob("*.Z.sac"))]

names = [header.station for header in headers]
coordinates = np.array([[header.sac.unused6, header.sac.unused7] for header in headers], dtype=float)

nstation = len(names)
dt       = headers[0].delta
npts     = headers[0].npts
nlag     = round(MAX_LAG / dt)
nwindow  = len(windows)
npair    = nstation * (nstation - 1) // 2

# 原始时间窗补零到至少 2*npts-1，使用同一 FFT 长度进行白化和互相关
nfft = next_fast_len(2 * npts - 1, real=True)

# 各时间窗的采样一致，白化频带渐变窗和时延索引只需确定一次
frequencies = np.fft.rfftfreq(nfft, d=dt)
band = (frequencies >= FREQMIN) & (frequencies <= FREQMAX)

taper = np.zeros_like(frequencies)
taper[band] = tukey(np.count_nonzero(band), alpha=2 * TAPER_FRACTION)

# 逆变换后负时延位于数组末尾，零时延和正时延位于数组开头
lag_indices = np.arange(-nlag, nlag + 1) % nfft


# 重建互相关目录
if OUTPUT.exists():
    shutil.rmtree(OUTPUT)

OUTPUT.mkdir()
for component in ("ZZ", "RR", "TT"):
    (OUTPUT / component).mkdir()


# 一次只处理一对台站，内存中仅保留当前时间窗和这一对的三个累加频谱
# 当前台站对的全部时间窗处理完后立即写出，再处理下一对
pair_index = 0

for source in range(nstation):
    for receiver in range(source + 1, nstation):
        pair_index += 1
        source_name   = names[source]
        receiver_name = names[receiver]

        # 方位由实际坐标计算，两个台站均以源台连线为径向
        north, east = coordinates[receiver] - coordinates[source]
        azimuth  = np.arctan2(east, north)
        distance = np.hypot(north, east)

        stacks = {component: np.zeros_like(taper, dtype=complex) for component in ("ZZ", "RR", "TT")}

        for window_index, window in enumerate(windows, 1):
            source_record   = read_segment(window, source_name)
            receiver_record = read_segment(window, receiver_name)

            # 正时延表示接收台站相对虚拟源延迟
            source_white   = process(source_record["Z"], nfft, taper)
            receiver_white = process(receiver_record["Z"], nfft, taper)

            stacks["ZZ"] += receiver_white * np.conj(source_white)

            # 水平分量先按当前台站对的方位旋转，再做与垂向相同的预处理
            radial_source, transverse_source = rotate_horizontal(source_record["N"], source_record["E"], azimuth)
            radial_receiver, transverse_receiver = rotate_horizontal(receiver_record["N"], receiver_record["E"], azimuth)

            for component, source_data, receiver_data in (
                    ("RR", radial_source, radial_receiver), ("TT", transverse_source, transverse_receiver)):
                source_white   = process(source_data, nfft, taper)
                receiver_white = process(receiver_data, nfft, taper)

                stacks[component] += receiver_white * np.conj(source_white)

            print(f"[{pair_index}/{npair}] {source_name} - {receiver_name}: window {window_index}/{nwindow}", end="\r", flush=True)

        # 频域逐窗等权平均，每个分量仅逆变换一次，再按负时延到正时延的顺序截取
        for component, stack in stacks.items():
            correlation = np.fft.irfft(stack / nwindow, n=nfft)
            result = Trace(correlation[lag_indices].astype(np.float32))
            result.stats.delta   = dt
            result.stats.channel = component
            result.stats.station = receiver_name

            # SAC 头段保存台站对和水平距离，供绘图脚本直接使用
            result.stats.sac = {
                "b":      -nlag * dt,
                "kevnm":  source_name,
                "dist":   float(distance),
            }

            result.write(str(OUTPUT / component / f"{source_name}_{receiver_name}.{component}.sac"), format="SAC")

        print(f"[{pair_index}/{npair}] {source_name} - {receiver_name}: saved ZZ, RR, TT ({nwindow} windows)", flush=True)
