import matplotlib.pyplot as plt
import pygrt

# BEGIN GEOMETRY
fig, axes = pygrt.utils.plot_sCoulomb_faults_geometry("fault.inp", outpath="fault.svg")
plt.close(fig)
# END GEOMETRY

# BEGIN FIELD PLOTS
data = pygrt.utils.read_nc_variables("stsyn_ff.nc")
fig, ax = pygrt.utils.plot_static_displacement(data, src_fault="fault.inp", clim_ratio=0.5, outpath="disp.svg")
plt.close(fig)
fig, axes = pygrt.utils.plot_static_coulomb(data, src_fault="fault.inp", clim_ratio=0.5, outpath="coulomb.svg")
plt.close(fig)

# END FIELD PLOTS
