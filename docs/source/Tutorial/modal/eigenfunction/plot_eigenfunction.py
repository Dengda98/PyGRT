import pygrt

egn_R = pygrt.utils.read_eigenfunction("egn_R.nc")
egn_L = pygrt.utils.read_eigenfunction("egn_L.nc")
pygrt.utils.plot_eigenfunction({"Rayleigh": egn_R, "Love": egn_L}, scale=2, outpath="eigenfunction.svg")
