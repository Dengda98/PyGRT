#!/bin/bash

set -euo pipefail

# BEGIN

M="milrow"
D="2/0"
R="100"
n="5"

# Calculate Full-wave
grt greenfn -M$M -D$D -R$R -N200/0.5+n$n -OGRN -e

# Calculate Surface-wave
grt eigenv -M$M -F0/1/0.01 -SR -N -Cphase_R.nc 
grt eigenv -M$M -F0/1/0.01 -SL -N -Cphase_L.nc 

# 0th order
grt modsum -Cphase_R.nc -D$D -R$R -N0 -OGRN_NM_0 -W$n -e
grt modsum -Cphase_L.nc -D$D -R$R -N0 -OGRN_NM_0 -W$n -e
# 1st order
grt modsum -Cphase_R.nc -D$D -R$R -N1 -OGRN_NM_1 -W$n -e
grt modsum -Cphase_L.nc -D$D -R$R -N1 -OGRN_NM_1 -W$n -e
# 2nd order
grt modsum -Cphase_R.nc -D$D -R$R -N2 -OGRN_NM_2 -W$n -e
grt modsum -Cphase_L.nc -D$D -R$R -N2 -OGRN_NM_2 -W$n -e
# all
grt modsum -Cphase_R.nc -D$D -R$R -N -OGRN_NM_all -W$n -e
grt modsum -Cphase_L.nc -D$D -R$R -N -OGRN_NM_all -W$n -e
# END

python run.py

for mode in 0 1 2 all; do
    python compare_wave.py "$mode"
done

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -rf GRN* *.nc
