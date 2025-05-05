import numpy as np
import matplotlib.pyplot as plt
import pygrt 


modarr = np.loadtxt("milrow")

pymod = pygrt.PyModel1D(modarr, 10, 0)

st_grn = pymod.compute_grn(
    distarr=[1800], 
    nt=1400, 
    dt=1, 
    Length=20, 
    safilonTol=1e-2,
    delayT0=100,
)[0]

st_grn.plot(equal_scale=False, outfile='test.pdf')