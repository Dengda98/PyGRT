"""对 ZZ 互相关做 FTAN，并与基阶 Rayleigh 理论频散对比"""

from pathlib import Path

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from scipy.fft import fft, fftfreq, ifft, next_fast_len
from scipy.signal import find_peaks
from scipy.signal.windows import tukey

import pygrt


# 在本脚本所在目录运行，使用 ex25 已有的互相关和理论频散文件
INPUT           = "ambient_ccf/ZZ/STA001_STA011.ZZ.sac"
PHASE_INPUT     = "ambient_run/phaseR.nc"
GROUP_INPUT     = "ambient_run/groupR.nc"
OUTPUT          = Path("ambient_figures")
FREQMIN         = 0.1                          # 预滤波频带下限，Hz
FREQMAX         = 1.5                          # 预滤波频带上限，Hz
FREQUENCIES     = np.linspace(FREQMIN, FREQMAX, 281)
VELOCITIES      = np.linspace(1.0, 4.5, 401)    # 群速度和相速度绘图网格，km/s
ALPHA           = 20.0                         # 群速度使用固定的高斯滤波参数
GAMMA           = np.sqrt(ALPHA / np.pi)           # 相速度使用随频率变化的高斯滤波宽度参数
GROUP_WINDOW    = (1.4, 3.5)                    # 群速度搜索范围，km/s
PHASE_REFERENCE = 3.2                          # 最低中心频率处的参考相速度，km/s，只用于确定追踪起点
PHASE_CYCLE     = 0                            # 相位整周数 n，本例取 0，实际数据需结合参考相速度等先验信息确定


def read_correlation(path):
    """带通滤波后将负时延支翻转，与正时延支等权平均得到单边信号

    :param path: 双边互相关的 SAC 文件
    """
    trace = read(path)[0]
    trace.data = trace.data.astype(float)
    trace.taper(max_percentage=0.05, type="cosine")
    trace.filter("bandpass", freqmin=FREQMIN, freqmax=FREQMAX, corners=4, zerophase=True)

    middle  = trace.stats.npts // 2
    stacked = (trace.data[middle:] + trace.data[middle::-1]) / 2

    times = np.arange(stacked.size) * trace.stats.delta
    return times, stacked, trace.stats.delta, float(trace.stats.sac.dist)


def compute_ftan(times, data, dt):
    """对同一条单边信号做高斯窄带滤波，返回包络和实波形

    :param times: 从零开始的单边时延采样，s
    :param data: 正负支翻转平均后的单边互相关
    :param dt: 采样间隔，s
    """
    nfft     = next_fast_len(4 * data.size)
    fft_freq = fftfreq(nfft, dt)
    window   = tukey(data.size, alpha=0.1)
    spectrum = 2 * fft(data * window, nfft)
    spectrum[fft_freq <= 0] = 0

    envelopes = np.empty((FREQUENCIES.size, data.size))
    waveforms = np.empty_like(envelopes)

    for index, f0 in enumerate(FREQUENCIES):
        band = spectrum * np.exp(-ALPHA * ((fft_freq - f0) / f0) ** 2)
        envelopes[index] = np.abs(ifft(band)[:data.size])

        alpha = GAMMA ** 2 * 2 * np.pi * f0
        band  = spectrum * np.exp(-alpha * ((fft_freq - f0) / f0) ** 2)
        waveforms[index] = ifft(band)[:data.size].real

    return times, envelopes, waveforms


def find_arrival_peaks(row, times):
    """寻找正的局部峰，并用三点抛物线插值细化峰到时

    :param row: 单个中心频率对应的包络或窄带实波形
    :param times: 单边时延采样，s
    """
    peaks = find_peaks(row)[0]
    peaks = peaks[row[peaks] > 0]
    left, center, right = row[peaks - 1], row[peaks], row[peaks + 1]
    offset   = (left - right) / (2 * (left - 2 * center + right))
    arrivals = times[peaks] + offset * (times[1] - times[0])
    return peaks, arrivals


def phase_time_shift(frequency, cycle):
    """将互相关的 pi/4 相位和整周数换算为到时修正，s

    :param frequency: 中心频率，Hz
    :param cycle: 相位整周数 n
    """
    return (1 / 8 + cycle) / frequency


def track_arrivals(values, times, distance, velocity_window, time_shift=0, start_time=None, reverse=False):
    """从低频向高频跟踪相邻峰，未指定起点时从速度窗内最强峰起步

    :param values: 按中心频率递增排列的包络或窄带实波形
    :param times: 正时延采样，s
    :param distance: 台间距，km
    :param velocity_window: 允许拾取的速度范围，km/s
    :param time_shift: 各中心频率对应的到时修正，s，群速度使用零修正
    :param start_time: 起始频率处的参考峰到时，s
    :param reverse: 是否改为从高频向低频追踪
    """
    arrivals   = np.full(FREQUENCIES.size, np.nan)
    time_shift = np.broadcast_to(time_shift, FREQUENCIES.shape)
    previous   = start_time
    indices    = range(FREQUENCIES.size)
    if reverse:
        indices = reversed(indices)

    for index in indices:
        row = values[index]
        peaks, refined = find_arrival_peaks(row, times)
        speeds = distance / (refined + time_shift[index])
        keep   = (speeds >= velocity_window[0]) & (speeds <= velocity_window[1])
        peaks, refined = peaks[keep], refined[keep]

        # 没有局部峰的频率不拾取
        if peaks.size == 0:
            continue

        selected = np.argmax(row[peaks]) if previous is None else np.argmin(np.abs(refined - previous))
        arrivals[index] = previous = refined[selected]

    return arrivals


def plot_ftan(times, envelopes, waveforms, distance, group_picks, phase_picks, phase_shift):
    """绘制群速度、相速度彩图，叠加拾取结果和基阶理论频散

    :param times: 正时延采样，s
    :param envelopes: 各中心频率对应的包络
    :param waveforms: 各中心频率对应的窄带实波形
    :param distance: 台间距，km
    :param group_picks: 群速度拾取，km/s
    :param phase_picks: 相速度拾取，km/s
    :param phase_shift: 起始定支后，各中心频率共用的相位到时修正，s
    """
    normalized_envelopes = envelopes / envelopes.max(axis=1, keepdims=True)
    normalized_waveforms = waveforms / np.max(np.abs(waveforms), axis=1, keepdims=True)
    group_image = np.array([np.interp(distance / VELOCITIES, times, row) for row in normalized_envelopes])
    phase_times = distance / VELOCITIES[None, :] - phase_shift[:, None]
    phase_image = np.array([np.interp(arrivals, times, row, left=np.nan, right=np.nan)
                            for arrivals, row in zip(phase_times, normalized_waveforms)])
    # 在显示速度窗内逐频率重新归一化相速度图
    phase_image /= np.nanmax(np.abs(phase_image), axis=1, keepdims=True)

    # 理论频散仅用于绘图对照，不参与拾取
    phase_reference = pygrt.io.read_dispersion(PHASE_INPUT)
    group_reference = pygrt.io.read_dispersion(GROUP_INPUT)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.3), layout="constrained")
    outline = [path_effects.withStroke(linewidth=3.5, foreground="white")]
    panels = (
        ("Group", group_image, group_picks, group_reference, "u", 0),
        ("Phase", phase_image, phase_picks, phase_reference, "c", -1),
    )
    for ax, panel in zip(axes, panels):
        kind, image, picked, reference, key, vmin = panel
        ax.pcolormesh(FREQUENCIES, VELOCITIES, image.T, cmap="jet", vmin=vmin, vmax=1, shading="auto", rasterized=True)
        ax.plot(FREQUENCIES, picked, color="#c93b32", linewidth=1.8, label="FTAN pick", path_effects=outline)
        ax.plot(reference["freq"][0], reference[key][0], color="#174b7a", linestyle="--", linewidth=1.8,
                label="Fundamental Rayleigh", path_effects=outline)
        ax.set(xlabel="Frequency (Hz)", ylabel=f"{kind} velocity (km/s)",
               xlim=(FREQMIN, FREQMAX), ylim=(VELOCITIES[0], VELOCITIES[-1]))
        ax.legend(loc="upper right", fontsize=8)

    fig.savefig(OUTPUT / "ftan.png", dpi=180)
    plt.close(fig)


OUTPUT.mkdir(exist_ok=True)

times, data, dt, distance = read_correlation(INPUT)
times, envelopes, waveforms = compute_ftan(times, data, dt)
group_arrivals = track_arrivals(envelopes, times, distance, GROUP_WINDOW, reverse=True)

# 将最低频率处的参考相速度换算为起始峰到时，再向高频追踪同一分支
phase_shift    = phase_time_shift(FREQUENCIES, PHASE_CYCLE)
phase_start    = distance / PHASE_REFERENCE - phase_shift[0]
phase_window   = (VELOCITIES[0], VELOCITIES[-1])
phase_arrivals = track_arrivals(waveforms, times, distance, phase_window, phase_shift, phase_start)

# 群速度使用中心频率，也可升级为群峰处的瞬时频率
group_picks = distance / group_arrivals
# 相速度用窄带实波形正峰到时加上统一的相位修正
phase_picks = distance / (phase_arrivals + phase_shift)
plot_ftan(times, envelopes, waveforms, distance, group_picks, phase_picks, phase_shift)

# DATA = Path("ambient_ftan")
# DATA.mkdir(exist_ok=True)
# np.savetxt(DATA / "ftan_picks.txt",
#            np.column_stack((FREQUENCIES, group_arrivals, group_picks, phase_arrivals, phase_picks)), fmt="%.8g",
#            header="f0_Hz group_arrival_s group_velocity_km_s phase_peak_s phase_velocity_km_s")
