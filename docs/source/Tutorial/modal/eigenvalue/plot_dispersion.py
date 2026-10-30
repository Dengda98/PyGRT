import pygrt

phase_R = pygrt.io.read_dispersion("phase_R.nc")
phase_L = pygrt.io.read_dispersion("phase_L.nc")
pygrt.plot.plot_dispersion({"Rayleigh": phase_R, "Love": phase_L}, outpath="dispersion_py.svg")
