from obspy import *
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8, 9))
subfigs = fig.subfigures(2, 1, hspace=0.1)

st = read("GRN/milrow_2_0_100/DS*")
st.plot(fig=subfigs[0])
for ax in subfigs[0].get_axes():
    ax.ticklabel_format(axis='y', style='sci', scilimits=[0,0])
subfigs[0].suptitle('No Upsampling')

st = read("GRN_up/milrow_2_0_100/DS*")
st.plot(fig=subfigs[1])
for ax in subfigs[1].get_axes():
    ax.ticklabel_format(axis='y', style='sci', scilimits=[0,0])
subfigs[1].suptitle('Upsampling = 5')

fig.savefig("upsampling.svg", bbox_inches='tight')