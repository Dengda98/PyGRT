from obspy import *
import matplotlib.pyplot as plt

def plot(pattern:str):
    st = read(pattern)
    fig = plt.figure(figsize=(8, len(st)))
    st.plot(fig=fig, equal_scale=False)
    for ax in fig.get_axes():
        ax.ticklabel_format(axis='y', style='sci', scilimits=[0,0])

    return fig

fig = plot("syn/?.sac")
fig.savefig("disp.svg", bbox_inches='tight')
fig = plot("syn/stress*.sac")
fig.savefig("stress.svg", bbox_inches='tight')