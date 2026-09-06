import pygrt

phase_R = pygrt.utils.read_dispersion("phase_R.nc")
phase_L = pygrt.utils.read_dispersion("phase_L.nc")
pygrt.utils.plot_dispersion({"Rayleigh": phase_R, "Love": phase_L}, outpath="dispersion_py.svg")
