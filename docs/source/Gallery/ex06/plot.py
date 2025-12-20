import numpy as np
import matplotlib.pyplot as plt
from obspy import *

st = read("./GRN_SAFIM/milrow_2_0_1800/*")

fig = plt.figure(figsize=(10, 0.8*len(st)))
st.plot(equal_scale=False, fig=fig, linewidth=0.5)

fig.savefig(f"safim.svg", bbox_inches='tight')