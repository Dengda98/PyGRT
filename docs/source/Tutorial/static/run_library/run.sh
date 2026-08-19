#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

rm -f stgrn.nc stsyn.nc

# BEGIN GRN
# -R 生成一维震中距列表，并与多深度组合写入一个 NetCDF 文件
grt static greenfn -Mmilrow -Ds2,4 -Dr0,2 -R0,5,10,15 -Ostgrn.nc
# END GRN

# BEGIN SYN
# 目标深度和二维接收网格不必是库中的直接采样值，static syn 会进行插值
grt static syn -Gstgrn.nc -Ds3 -Dr1 -S1e24 -M33/50/120 -X-10/10/5 -Y-10/10/5 -Ostsyn.nc
# END SYN

rm -f stgrn.nc stsyn.nc
