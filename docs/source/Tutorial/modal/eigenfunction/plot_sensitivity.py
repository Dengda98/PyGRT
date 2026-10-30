import matplotlib.pyplot as plt
import pygrt

for char in ["c", "u"]:
    for wtype, title in [("R", "Rayleigh"), ("L", "Love")]:
        data = pygrt.io.read_sensitivity(f"{char}sens_{wtype}.nc")
        fig, _ = pygrt.plot.plot_sensitivity(data, title=title, outpath=f"{char}sens_{wtype}.svg")
        plt.close(fig)
