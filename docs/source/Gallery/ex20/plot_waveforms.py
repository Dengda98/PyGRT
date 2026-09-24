import math

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pygrt
from scipy.optimize import brentq


NU = 0.25
R = 20.0
SOURCE_DEPTH = 5.0
RECEIVER_DEPTH = 4.0
AZIMUTH = 0.0
DTBAR = 0.002
TMAX = 2.0
COMPONENT = (2, 2)
PHASES = ("P", "S", "PP", "SS", "PS", "SP", "sPs")
K = math.sqrt(0.5 * (1.0 - 2.0 * NU) / (1.0 - NU))
THETA_C = math.asin(K)


def conversion_split(p_depth, s_depth):
    """Return the surface conversion point for a P-S path."""

    def residual(horizontal_split):
        p_path = math.hypot(horizontal_split, p_depth)
        s_path = math.hypot(R - horizontal_split, s_depth)
        return (R - horizontal_split) / s_path - K * horizontal_split / p_path

    left = 1e-12 * (1.0 + R)
    right = R - left
    return brentq(residual, left, right)


def phase_arrivals():
    direct_distance = math.hypot(R, SOURCE_DEPTH - RECEIVER_DEPTH)
    reflected_distance = math.hypot(R, SOURCE_DEPTH + RECEIVER_DEPTH)
    theta_reflection = math.atan2(R, SOURCE_DEPTH + RECEIVER_DEPTH)
    ps_split = conversion_split(SOURCE_DEPTH, RECEIVER_DEPTH)
    sp_split = conversion_split(RECEIVER_DEPTH, SOURCE_DEPTH)
    return {
        "P": K,
        "S": 1.0,
        "PP": K * reflected_distance / direct_distance,
        "SS": reflected_distance / direct_distance,
        "PS": (K * math.hypot(ps_split, SOURCE_DEPTH) + math.hypot(R - ps_split, RECEIVER_DEPTH)) / direct_distance,
        "SP": (K * math.hypot(sp_split, RECEIVER_DEPTH) + math.hypot(R - sp_split, SOURCE_DEPTH)) / direct_distance,
        "sPs": math.cos(theta_reflection - THETA_C) * reflected_distance / direct_distance,
    }


tbar = np.arange(0.0, TMAX + 0.5 * DTBAR, DTBAR)
phase_data = {}
for phase in PHASES:
    green, _, _, _ = pygrt.utils.lamb3(
        nu=NU,
        tbar=tbar,
        R=R,
        depsrc=SOURCE_DEPTH,
        deprcv=RECEIVER_DEPTH,
        azimuth=AZIMUTH,
        phases=[phase],
    )
    phase_data[phase] = green[:, COMPONENT[0], COMPONENT[1]]

sum_data = sum(phase_data.values())
arrivals = phase_arrivals()
amplitude = max(np.max(np.abs(data)) for data in (*phase_data.values(), sum_data))
amplitude = max(amplitude, 1e-12) * 0.3

trace_names = (*PHASES, "Sum")
trace_data = (*[phase_data[phase] for phase in PHASES], sum_data)
offsets = np.arange(len(trace_names) - 1, -1, -1, dtype=float) * 1.6

fig, axis = plt.subplots(figsize=(7, 8))
for offset, phase, data in zip(offsets[:-1], PHASES, trace_data[:-1]):
    axis.plot(tbar, data / amplitude + offset, color="#1f4e79", linewidth=1.5, clip_on=False)
    axis.text(tbar[0]-0.02, offset, phase, color="black", fontsize=13, ha="right", va="center", clip_on=False)
    arrival = arrivals[phase]
    axis.vlines(arrival, offset - 0.52, offset + 0.52, color="#c0392b", linestyle="--", linewidth=0.8)
    axis.text(
        arrival + 0.012,
        offset + 0.55,
        f"{phase}  {arrival:.3f}",
        color="#c0392b",
        fontsize=9,
        ha="left",
        va="bottom",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5),
    )

axis.plot(tbar, sum_data / amplitude + offsets[-1], color="black", linewidth=2.0, clip_on=False)
axis.text(tbar[0]-0.02, offsets[-1], "Sum", color="black", fontsize=13, ha="right", va="center", clip_on=False)
axis.set_xlim(tbar[0], tbar[-1])
axis.set_ylim(-0.7, offsets[0] + 0.7)
axis.set_yticks([])
for spine in axis.spines.values():
    spine.set_visible(False)
axis.spines["bottom"].set_visible(True)
axis.tick_params(axis="x", bottom=True, labelbottom=True)
axis.set_xlabel(r"Normalized time $\bar{t}=\beta t/r$")
axis.set_title(r"Lamb3 phase-separated waveforms: $G_{33}$", pad=15)
fig.subplots_adjust(left=0.08, right=0.99, bottom=0.10, top=0.93)
fig.savefig("phase_waveforms.svg", bbox_inches="tight")
plt.close(fig)
