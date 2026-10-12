#!/bin/bash
set -euo pipefail

# BEGIN DISPERSION
# Calculate Rayleigh-wave dispersion
grt eigenv -Mmod2 -SR -F0/5/0.01 -N -Cphase_R.nc
# Calculate Love-wave dispersion
grt eigenv -Mmod2 -SL -F0/5/0.01 -N -Cphase_L.nc
# END DISPERSION

# BEGIN EIGENFN
# Calculate Rayleigh-wave eigenfunction
grt eigenfn -Cphase_R.nc -F1 -N -Wegn_R.nc+z0/50/0.1
# Calculate Love-wave eigenfunction
grt eigenfn -Cphase_L.nc -F1 -N -Wegn_L.nc+z0/50/0.1
# END EIGENFN

# BEGIN GROUP
# Calculate Rayleigh-wave group velocity
grt eigenfn -Cphase_R.nc -N -Ugroup_R.nc
# Calculate Love-wave group velocity
grt eigenfn -Cphase_L.nc -N -Ugroup_L.nc
# END GROUP

# BEGIN SENS
# Calculate Rayleigh-wave dispersion sensitivity
grt eigenfn -Cphase_R.nc -F0/0.5/0.1 -N0 -K+ccsens_R.nc+uusens_R.nc+z0.2
# Calculate Love-wave dispersion sensitivity
grt eigenfn -Cphase_L.nc -F0/0.5/0.1 -N0 -K+ccsens_L.nc+uusens_L.nc+z0.2
# END SENS

python run.py
python plot_eigenfunction.py
python plot_dispersion_group.py
python plot_sensitivity.py

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -f *.nc
