import pygrt
import numpy as np
import matplotlib.pyplot as plt


data1, data2, ptamdata, dist = pygrt.utils.read_statsfile_ptam("GRN_grtstats/halfspace2_0.01_0/PTAM_0000_1.00000e+00/PTAM_0142_*")
fig, ax = pygrt.utils.plot_statsdata_ptam(data1, data2, ptamdata, dist, "DC", "2", "2")
# fig.tight_layout()
# fig.savefig("view_stats.png", dpi=100, bbox_inches='tight')
# plt.show()
