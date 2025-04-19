import pygrt
import numpy as np
import matplotlib.pyplot as plt
import glob

data = pygrt.utils.read_statsfile("stgrtstats/halfspace2_0.1_0/K")

dist = 3
fig, ax = pygrt.utils.plot_statsdata(data, dist, "DC", "2", "2")


data1, data2, ptamdata, dist = pygrt.utils.read_statsfile_ptam("stgrtstats/halfspace2_0.1_0/PTAM_0090_*/PTAM")
fig, ax = pygrt.utils.plot_statsdata_ptam(data1, data2, ptamdata, dist, "DC", "2", "2")
plt.show()
