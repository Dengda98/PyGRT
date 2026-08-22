#!/bin/bash

set -euo pipefail

rm -f stgrn.nc stsyn.nc

# BEGIN GRN
# -R 生成一维震中距列表，并与多深度组合写入一个 NetCDF 文件
# 建议使用 -e，这样计算的格林函数库包含位移空间偏导数，可以用于后续相关的计算
grt static greenfn -Mmilrow -Ds2,4 -Dr0,2 -R0,5,10,15 -Ostgrn.nc -e
# END GRN

# BEGIN SYN
# 使用 -Ds, -Dr 明确给出点源的震源深度和台站深度
# 使用 -e 表示增加位移偏导数的合成
grt static syn -Gstgrn.nc -Ds3 -Dr1 -S1e24 -M33/50/120 -X-10/10/5 -Y-10/10/5 -Ostsyn.nc -e
# END SYN

rm -f stgrn.nc stsyn.nc
