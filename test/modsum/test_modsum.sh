#!/bin/bash

set -euo pipefail

# 计算面波频散, eigenv -F 和 greenfn -N 相匹配
grt eigenv -M../milrow -F0/1/0.01 -SR -N -Cphase_R.nc
grt eigenv -M../milrow -F0/1/0.01 -SL -N -Cphase_L.nc

# 模态叠加得到格林函数
# 仅 0 阶
grt modsum -Cphase_R.nc -D2/1 -R100 -N0 -OGRN_NM_0 -W5 -e
grt modsum -Cphase_L.nc -D2/1 -R100 -N0 -OGRN_NM_0 -W5 -e
# 仅 1 阶
grt modsum -Cphase_R.nc -D2/1 -R100 -N1 -OGRN_NM_1 -W5 -e
grt modsum -Cphase_L.nc -D2/1 -R100 -N1 -OGRN_NM_1 -W5 -e
# 仅 2 阶
grt modsum -Cphase_R.nc -D2/1 -R100 -N2 -OGRN_NM_2 -W5 -e
grt modsum -Cphase_L.nc -D2/1 -R100 -N2 -OGRN_NM_2 -W5 -e
# 全部
grt modsum -Cphase_R.nc -D2/1 -R100 -N -OGRN_NM_all -W5 -e
grt modsum -Cphase_L.nc -D2/1 -R100 -N -OGRN_NM_all -W5 -e

rm -rf GRN* *.nc