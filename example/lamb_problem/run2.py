import numpy as np
from obspy import *
import matplotlib.pyplot as plt
import pygrt

# load halfspace model
modarr = np.loadtxt("halfspace")
_, Vp, Vs, Rho, _, _ = modarr
modarr = modarr[None,:]

# shear module
mu = Vs**2 * Rho

depsrc = 0.0
deprcv = 0.0

rs = np.array([10])

nt = 500
dt = 0.01 

r = rs[0]

t = np.arange(0, nt)*dt * Vs/r 

pymod = pygrt.PyModel1D(modarr, depsrc, deprcv)

# compute Green's Functions
st_grn = pymod.gen_gf_spectra(
    distarr=rs, 
    nt=nt, 
    dt=dt, 
    keps=-1,
    Length=20, 
    gf_source=['VF', 'HF']
)[0]

# convolve with step function
pygrt.utils.stream_integral(st_grn)

st_grn.plot(equal_scale=False, figsize=(8, 6), outfile="GRN2.pdf", show=False)
