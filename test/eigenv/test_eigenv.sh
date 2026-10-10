#!/bin/bash

set -euo pipefail
source ../common.sh

grt eigenv -h
grt disp2asc -h
grt eigenv -M../milrow -SR -F0/0.5/0.1 -N -Cphase_R.nc
grt eigenv -M../milrow -SL -F0/0.5/0.1 -N2 -T+t3 -P2 -Cphase_L.nc
grt eigenv -M../milrow -SR -F2/4/1+p -Cphase_p.nc
grt eigenv -M../milrow -SR -X1+c3/5.5+i0 > secular.txt
grt disp2asc -Cphase_R.nc -N > phase.txt
grt disp2asc -Cphase_R.nc -F0.2/0.5/0.1 -N0/1 > phase.txt
grt disp2asc -Cphase_p.nc -F2+p > phase.txt

# 面波类型不可设置为 R、L 之外的值
expect_fail grt eigenv -M../milrow -SX -F0.5 -Cphase_bad.nc
# 频率采样间隔不可为零或负数
expect_fail grt eigenv -M../milrow -SR -F0/1/0 -Cphase_bad.nc
# 不可同时设置相速度文件和群速度文件
expect_fail grt disp2asc -Cphase_R.nc -Uphase_L.nc

# 周期不可为零或负数
expect_fail grt eigenv -M../milrow -SR -F0+p -Cphase_bad.nc

python -u test_eigenv.py

rm -rf phase_*.nc phase.txt secular.txt
