import pygrt

sec_R = pygrt.utils.read_secfunc("secfunc_R")
sec_L = pygrt.utils.read_secfunc("secfunc_L")
pygrt.utils.plot_secfunc({"Rayleigh": sec_R[0], "Love": sec_L[0]}, outpath="secfunc.svg")
