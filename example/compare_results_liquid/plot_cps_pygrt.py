import pygrt
import numpy as np 
import matplotlib.pyplot as plt 
from obspy import * 


# load model 
modarr = np.loadtxt("./liquid")

depsrc = 20.0   
deprcv = 8.2    
pymod = pygrt.PyModel1D(modarr, depsrc, deprcv) 

rs = np.array([100]) 

nt = 512   
dt = 0.2   
zeta = 0.8  

# compute Green's Functions
st_grt = pymod.compute_grn(
    distarr=rs, 
    nt=nt, 
    dt=dt, 
    zeta=zeta, 
    Length=20,
)[0]


try:
    st_cps = read("liquid_sdep20_rdep8.2/GRN/*")

    from plot_cps_grt import plot
    plot(st_grt, st_cps, "compare_cps_pygrt.pdf")
except Exception as e:
    print(str(e))
