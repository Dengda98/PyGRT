import pygrt

group_R = pygrt.io.read_dispersion("group_R.nc")
group_L = pygrt.io.read_dispersion("group_L.nc")
pygrt.plot.plot_dispersion({"Rayleigh": group_R, "Love": group_L}, outpath="dispersion_py_group.svg")
