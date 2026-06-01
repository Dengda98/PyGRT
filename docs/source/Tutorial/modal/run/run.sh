#!/bin/bash

# BEGIN DISPERSION
# Calculate Rayleigh-wave dispersion
grt eigenv -Mmod2 -SR -F0/5/0.01 -N -Cphase_R.nc
# Calculate Love-wave dispersion
grt eigenv -Mmod2 -SL -F0/5/0.01 -N -Cphase_L.nc
# END DISPERSION

python plot_dispersion.py
bash plot_dispersion.sh

# BEGIN SECFUNC
freq=1
# Output Rayleigh-wave secular function
grt eigenv -Mmod2 -SR -X$freq+c3.0/5.5+i0  > secfunc_R
# Output Love-wave secular function
grt eigenv -Mmod2 -SL -X$freq+c3.0/5.5+i0  > secfunc_L
# END SECFUNC

python plot_secfunc.py

# BEGIN EIGENFN
# Calculate Rayleigh-wave eigenfunction
grt eigenfn -Cphase_R.nc -F1 -N -Wegn_R.nc+z0/50/0.1
# Calculate Love-wave eigenfunction
grt eigenfn -Cphase_L.nc -F1 -N -Wegn_L.nc+z0/50/0.1
# END EIGENFN

python plot_eigenfunction.py

# BEGIN GROUP
# Calculate Rayleigh-wave group velocity
grt eigenfn -Cphase_R.nc -N -Ugroup_R.nc
# Calculate Love-wave group velocity
grt eigenfn -Cphase_L.nc -N -Ugroup_L.nc
# END GROUP

python plot_dispersion_group.py
bash   plot_dispersion_group.sh

# BEGIN SENS
# Calculate Rayleigh-wave dispersion sensitivity
grt eigenfn -Cphase_R.nc -F0/0.5/0.1 -N0 -K+ccsens_R.nc+uusens_R.nc+z0.2
# Calculate Love-wave dispersion sensitivity
grt eigenfn -Cphase_L.nc -F0/0.5/0.1 -N0 -K+ccsens_L.nc+uusens_L.nc+z0.2
# END SENS

python plot_sensitivity.py


rm *.nc secfunc_*