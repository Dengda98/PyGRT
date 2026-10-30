import pygrt

sec_R = pygrt.io.read_secfunc("secfunc_R")
sec_L = pygrt.io.read_secfunc("secfunc_L")
pygrt.plot.plot_secfunc({"Rayleigh": sec_R[0], "Love": sec_L[0]}, outpath="secfunc.svg")
