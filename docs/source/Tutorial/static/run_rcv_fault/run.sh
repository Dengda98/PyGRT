#!/bin/bash

set -euo pipefail

rm -rf *.nc *.svg

# ------------------------------------------------
# BEGIN GRN
# 设置 -Ds0 仅考虑 5 km 的震源深度
# 设置 -Dr 考虑多个接收点深度
# 设置 -R 考虑多个震中距
# 设置 -e 增加位移偏导数的计算
grt static greenfn -Mmilrow -Ds5 -Dr0/10/0.5 -R0/20/0.5 -e -Ostgrn.nc
# END GRN
# ------------------------------------------------

# ------------------------------------------------
# BEGIN SYN
# 使用 -R 读取 Coulomb 格式的有限断层指定接收点
grt static syn -Gstgrn.nc -S1e24 -M33/90/0 -Rrcv_fault.inp -N -Ostsyn_rf.nc -e
# END SYN
# ------------------------------------------------

python plot_rcv_fault.py

rm -rf *.nc
