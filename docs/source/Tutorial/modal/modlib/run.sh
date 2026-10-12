#!/bin/bash
set -euo pipefail

rm -rf GRN phase_R.nc phase_L.nc

# BEGIN LIBRARY
# eigenv 的 -F 使用等间隔频率，供 modsum 进行逆傅里叶变换
grt eigenv -Mmilrow -F0/1/0.02 -SR -N -Cphase_R.nc
grt eigenv -Mmilrow -F0/1/0.02 -SL -N -Cphase_L.nc

# 计算全部震源深度、台站深度和震中距的组合
grt modsum -Cphase_R.nc -Ds2,4 -Dr0,2 -R80,100,120 -N0 -OGRN -W4
grt modsum -Cphase_L.nc -Ds2,4 -Dr0,2 -R80,100,120 -N0 -OGRN -W4
# END LIBRARY
