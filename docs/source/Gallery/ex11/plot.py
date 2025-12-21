from obspy import *
import matplotlib.pyplot as plt

def plot(pattern:str):
    st = read(pattern)
    fig = plt.figure(figsize=(6, 2*len(st)))
    st.plot(fig=fig, linewidth=0.5)
    for ax in fig.get_axes():
        ax.ticklabel_format(axis='y', style='sci', scilimits=[0,0])

    return fig

fig = plot("syn/[ZR].sac")
fig.savefig("disp.svg", bbox_inches='tight')