#!/bin/bash

set -euo pipefail
source ../common.sh

grt modsum -h
grt eigenv -M../milrow -SR -F0/0.5/0.1 -N -Cphase_R.nc
grt eigenv -M../milrow -SL -F0/0.5/0.1 -N -Cphase_L.nc
grt modsum -Cphase_R.nc -D2/1 -R100 -N0 -W2 -e -OGRN
grt modsum -Cphase_L.nc -D2/1 -R100 -N -W2 -e -OGRN
grt modsum -Cphase_R.nc -Ds1,2 -Dr0,1 -R100,120 -N0 -F0.1/0.4 -E-2/9 -P2 -OGRN
grt modsum -Cphase_R.nc -D2/1 -R100 -N0 -Ep-1 -Gev -OGRN

# 未指定多深度参数对应的接收深度
expect_fail grt modsum -Cphase_R.nc -Ds1,2 -R100 -N0 -OGRN
# 震中距列表不可包含逆序或重复值
expect_fail grt modsum -Cphase_R.nc -D2/1 -R100,50 -N0 -OGRN

# 升采样倍数不可为零或负数
expect_fail grt modsum -Cphase_R.nc -D2/1 -R100 -N0 -W0 -OGRN

python -u test_modsum.py

rm -rf GRN phase_R.nc phase_L.nc
