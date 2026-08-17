#!/usr/bin/env python3

from pathlib import Path
from typing import Dict, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.io import netcdf_file


HERE = Path(__file__).resolve().parent
CHANNELS = ("Z", "N", "E")


def read_grid(path: Path) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    """读取规则网格中的坐标和位移变量"""

    with netcdf_file(path, mode="r", mmap=False) as dataset:
        north = np.array(dataset.variables["north"][:], copy=True)
        east = np.array(dataset.variables["east"][:], copy=True)
        data = {name: np.array(dataset.variables[name][:], copy=True) for name in CHANNELS}
    return north, east, data


def plot_displacement(path: Path, output: Path, title: str) -> None:
    """绘制三个位移分量"""

    north, east, data = read_grid(path)
    extent = (east.min(), east.max(), north.min(), north.max())
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.4), layout='constrained')
    for ax, channel in zip(axes, CHANNELS):
        field = data[channel]
        limit = np.max(np.abs(field))
        image = ax.imshow(field, extent=extent, origin="lower", aspect="equal", cmap="seismic", vmin=-limit, vmax=limit)
        ax.set_title(channel)
        ax.set_xlabel("East (km)")
        ax.set_ylabel("North (km)")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="3.5%", pad=0.45)
        colorbar = fig.colorbar(image, cax=cax, orientation="horizontal")
        colorbar.set_label("Displacement (cm)", fontsize=8, labelpad=2)
        colorbar.ax.tick_params(labelsize=7, pad=1)
    fig.suptitle(title, y=1.05)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


plot_displacement(HERE / "okada_point.nc", HERE / "okada_point.svg", "Okada point source")
plot_displacement(HERE / "okada_finite.nc", HERE / "okada_finite.svg", "Okada finite fault")
