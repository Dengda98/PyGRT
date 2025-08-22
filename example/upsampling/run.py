import pygrt
import numpy as np 
import matplotlib.pyplot as plt 
from obspy import * 


# load model 
modarr = np.loadtxt("./milrow")

depsrc = 1.0   
deprcv = 0     
pymod = pygrt.PyModel1D(modarr, depsrc, deprcv) 

rs = np.array([100]) 

nt = 200   
dt = 0.5   

# compute Green's Functions
# uncomment "upsampling_n" to use upsampling
st_grt = pymod.compute_grn(
    distarr=rs, 
    nt=nt, 
    dt=dt, 
    # upsampling_n=5
)[0]

# fig = plt.figure(figsize=(10, 10))
# st_grt.plot(fig=fig, equal_scale=False)

# plt.show()