import pygrt

group_R = pygrt.utils.read_dispersion("group_R.nc")
group_L = pygrt.utils.read_dispersion("group_L.nc")
pygrt.utils.plot_dispersion({"Rayleigh": group_R, "Love": group_L}, outpath="dispersion_py_group.svg")
