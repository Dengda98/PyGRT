#!/bin/bash

set -euo pipefail
source ../common.sh

grt eigenfn -h
grt eigenv -M../milrow -SR -F0/0.5/0.1 -N -Cphase.nc
grt eigenfn -Cphase.nc -F0.5 -N0/2 -Wegn.nc+z0/10/1 -Ugroup.nc -K+ccsens.nc+uusens.nc+z1+xenergy.nc
grt eigenfn -Cphase.nc -F2+p -N0 -Wegn.nc+z1
grt eigenfn -Cphase.nc -F0.2/0.5 -N -Ugroup.nc
grt disp2asc -Ugroup.nc -N > group.txt

# 深度采样间隔不可为零或负数
expect_fail grt eigenfn -Cphase.nc -Wegn.nc+z0/10/0
# 输入的频散文件不存在
expect_fail grt eigenfn -Cmissing.nc -Ugroup.nc
# 频率采样间隔不可为零或负数
expect_fail grt disp2asc -Ugroup.nc -F0.1/0.5/0

# 模态阶数不可为负数
expect_fail grt eigenfn -Cphase.nc -N-1 -Ugroup.nc

python -u test_eigenfn.py

rm -rf phase.nc egn.nc group.nc csens.nc usens.nc energy.nc group.txt
