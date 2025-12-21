import numpy as np
import matplotlib.pyplot as plt
from obspy import *
import glob, os

def plot_one(prefix:str):
    st = read(f"./syn_dc/{prefix}*")
    for tr in st:
        if np.all(tr.data == 0.0):
            tr.data[:] = 1e-30
    fig = plt.figure(figsize=(8, len(st)))
    st.plot(fig=fig)
    for ax in fig.get_axes():
        ax.ticklabel_format(axis='y', style='sci', scilimits=[0,0])
    fig.savefig(f"{prefix}.svg", bbox_inches='tight')

plot_one("stress")
plot_one("strain")
plot_one("rotation")