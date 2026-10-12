#!/bin/bash
set -euo pipefail

# BEGIN DISPERSION
# Calculate Rayleigh-wave dispersion
grt eigenv -Mmod2 -SR -F0/5/0.01 -N -Cphase_R.nc
# Calculate Love-wave dispersion
grt eigenv -Mmod2 -SL -F0/5/0.01 -N -Cphase_L.nc
# END DISPERSION

# BEGIN SECFUNC
freq=1
# Output Rayleigh-wave secular function
grt eigenv -Mmod2 -SR -X$freq+c3.0/5.5+i0  > secfunc_R
# Output Love-wave secular function
grt eigenv -Mmod2 -SL -X$freq+c3.0/5.5+i0  > secfunc_L
# END SECFUNC

python run.py
python plot_dispersion.py
python plot_secfunc.py

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -f *.nc secfunc_*
