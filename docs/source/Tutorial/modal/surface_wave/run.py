import pygrt

# BEGIN
pymod = pygrt.PyModel1D(grn="GRN", modelpath="milrow")

# Calculate Full-wave
pymod.greenfn(depsrc=2.0, deprcv=0.0, dists=100.0, nt=200, dt=0.5, upsampling_n=5, calc_upar=True)

# Calculate Surface-wave
pymod.eigenv(wtype="R", freqs=(0.0, 1.0, 0.01), phase_path="phase_R.nc", all_modes=True)
pymod.eigenv(wtype="L", freqs=(0.0, 1.0, 0.01), phase_path="phase_L.nc", all_modes=True)

# 0th order
pymod.modsum(phase_path="phase_R.nc", depsrc=2.0, deprcv=0.0, dists=100.0, modes=0, output_path="GRN_NM_0", upsampling_n=5, calc_upar=True)
pymod.modsum(phase_path="phase_L.nc", depsrc=2.0, deprcv=0.0, dists=100.0, modes=0, output_path="GRN_NM_0", upsampling_n=5, calc_upar=True)
# 1st order
pymod.modsum(phase_path="phase_R.nc", depsrc=2.0, deprcv=0.0, dists=100.0, modes=1, output_path="GRN_NM_1", upsampling_n=5, calc_upar=True)
pymod.modsum(phase_path="phase_L.nc", depsrc=2.0, deprcv=0.0, dists=100.0, modes=1, output_path="GRN_NM_1", upsampling_n=5, calc_upar=True)
# 2nd order
pymod.modsum(phase_path="phase_R.nc", depsrc=2.0, deprcv=0.0, dists=100.0, modes=2, output_path="GRN_NM_2", upsampling_n=5, calc_upar=True)
pymod.modsum(phase_path="phase_L.nc", depsrc=2.0, deprcv=0.0, dists=100.0, modes=2, output_path="GRN_NM_2", upsampling_n=5, calc_upar=True)
# all
pymod.modsum(phase_path="phase_R.nc", depsrc=2.0, deprcv=0.0, dists=100.0, all_modes=True, output_path="GRN_NM_all", upsampling_n=5, calc_upar=True)
pymod.modsum(phase_path="phase_L.nc", depsrc=2.0, deprcv=0.0, dists=100.0, all_modes=True, output_path="GRN_NM_all", upsampling_n=5, calc_upar=True)
# END
