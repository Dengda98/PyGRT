import pygrt

# BEGIN DISPERSION
pymod = pygrt.PyModel1D(modelpath="mod2")
pymod.eigenv(wtype="R", freqs=(0.0, 5.0, 0.01), phase_path="phase_R.nc", all_modes=True)
pymod.eigenv(wtype="L", freqs=(0.0, 5.0, 0.01), phase_path="phase_L.nc", all_modes=True)
# END DISPERSION


# BEGIN SECFUNC
from contextlib import redirect_stdout

# 将输出重定向到文件中
with open("secfunc_R", "w") as f, redirect_stdout(f):
    pymod.eigenv(wtype="R", secular_freq=1.0, cmin=3.0, cmax=5.5, iref=0)
with open("secfunc_L", "w") as f, redirect_stdout(f):
    pymod.eigenv(wtype="L", secular_freq=1.0, cmin=3.0, cmax=5.5, iref=0)
# END SECFUNC


# BEGIN EIGENFN
pymod.eigenfn(
    phase_path="phase_R.nc",
    freqs=1.0,
    all_modes=True,
    eigenfn_path="egn_R.nc",
    depths=(0.0, 50.0, 0.1),
)
pymod.eigenfn(
    phase_path="phase_L.nc",
    freqs=1.0,
    all_modes=True,
    eigenfn_path="egn_L.nc",
    depths=(0.0, 50.0, 0.1),
)
# END EIGENFN


# BEGIN GROUP
pymod.eigenfn(phase_path="phase_R.nc", all_modes=True, group_path="group_R.nc")
pymod.eigenfn(phase_path="phase_L.nc", all_modes=True, group_path="group_L.nc")
# END GROUP


# BEGIN SENS
pymod.eigenfn(
    phase_path="phase_R.nc",
    freqs=(0.0, 0.5, 0.1),
    modes=0,
    csens_path="csens_R.nc",
    usens_path="usens_R.nc",
    sensitivity_dz=0.2,
)
pymod.eigenfn(
    phase_path="phase_L.nc",
    freqs=(0.0, 0.5, 0.1),
    modes=0,
    csens_path="csens_L.nc",
    usens_path="usens_L.nc",
    sensitivity_dz=0.2,
)
# END SENS
