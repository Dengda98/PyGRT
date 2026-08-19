#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

rm -rf GRN *.nc syn

# BEGIN LIBRARY
# eigenv 的 -F 使用等间隔频率，供 modsum 进行逆傅里叶变换
grt eigenv -Mmilrow -F0/1/0.02 -SR -N -Cphase_R.nc
grt eigenv -Mmilrow -F0/1/0.02 -SL -N -Cphase_L.nc

# -Ds/-Dr 和 -R 传入列表，modsum 在内部遍历全部组合
# 加上 -e 表示计算位移格林函数的空间偏导
grt modsum -Cphase_R.nc -Ds2,4 -Dr0,2 -R80,100,120 -N0 -OGRN -W4 -e
grt modsum -Cphase_L.nc -Ds2,4 -Dr0,2 -R80,100,120 -N0 -OGRN -W4 -e
# END LIBRARY

# BEGIN SYN
# syn 从面波格林函数库根目录精确选择一个深度和距离
# 加上 -e 表示计算合成位移的空间偏导
grt syn -GGRN -Ds4 -Dr2 -R100 -S1e24 -A30 -M33/50/120 -Osyn -e
# END SYN

rm -rf GRN *.nc syn
