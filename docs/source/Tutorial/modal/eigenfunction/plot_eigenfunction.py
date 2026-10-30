import pygrt

egn_R = pygrt.io.read_eigenfunction("egn_R.nc")
egn_L = pygrt.io.read_eigenfunction("egn_L.nc")
pygrt.plot.plot_eigenfunction({"Rayleigh": egn_R, "Love": egn_L}, scale=2, outpath="eigenfunction.svg")
