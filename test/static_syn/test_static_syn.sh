#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt static syn -h
grt static greenfn -M../milrow -D2/0 -R0/12/3 -e -Ostgrn.nc

# 点源机制、源强、坐标分量和位移导数
grt static syn -Gstgrn.nc -S1e20 -Ostsyn.nc
grt static syn -Gstgrn.nc -S1e20 -F2/-1/4 -Ostsyn.nc
grt static syn -Gstgrn.nc -S1e20 -M33/44/55 -Ostsyn.nc
grt static syn -Gstgrn.nc -Su1e6 -M33/44 -Ostsyn.nc
grt static syn -Gstgrn.nc -S1e20 -T1/-2/-5/0.5/3/1.2 -N -e -Ostsyn.nc
grt static syn -Gstgrn.nc -S1e20 -X-2/2/2 -Y-2/2/2 -Ostsyn.nc

# 多深度库、任意接收点、有限震源和接收断层
grt static greenfn -M../milrow -Ds1,2 -Dr0,1 -R0/12/3 -e -Ostgrn_multi.nc
grt static syn -Gstgrn_multi.nc -S1e20 -Ds1.5 -Dr0.5 -X-2/2/2 -Y-2/2/2 -P2 -Ostsyn.nc
grt static syn -Gstgrn_multi.nc -S1e20 -Ds2 -Qrcv_points.txt -Ostsyn.nc
grt static syn -Gstgrn_multi.nc -S1e20 -Ds2 -Qrcv_geometry.txt -N -e -Ostsyn.nc
grt static syn -Gstgrn_multi.nc -S1e20 -Ds2 -Urcv_faults.inr -Ostsyn.nc
grt static syn -Gstgrn_multi.nc -S1e20 -Ds2 -Urcv_faults.inr+i1/1 -Ostsyn.nc
grt static syn -Gstgrn_multi.nc -Cfaults.inp+i1/1 -Dr0 -X-2/2/2 -Y-2/2/2 -e -Ostsyn.nc
grt static syn -Gstgrn_multi.nc -Cfaults.inr -Qrcv_points.txt -Ostsyn.nc

# 显式不再剖分，每条矩形源及接收断层只取一个中心点
grt static syn -Gstgrn_multi.nc -Cfaults.inp+i0/0 -Urcv_faults.inr+i0/0 -e -Ostsyn.nc

# 未指定多深度库对应的震源深度和接收深度
expect_fail grt static syn -Gstgrn_multi.nc -S1e20 -Ostsyn.nc
# 震源深度不可超出格林函数库范围
expect_fail grt static syn -Gstgrn_multi.nc -S1e20 -Ds3 -Dr0 -Ostsyn.nc
# 不可同时设置有限震源文件和点源源强
expect_fail grt static syn -Gstgrn_multi.nc -S1e20 -Cfaults.inp -Ostsyn.nc
# 不可同时设置接收点文件和接收网格
expect_fail grt static syn -Gstgrn_multi.nc -S1e20 -Ds2 -Qrcv_points.txt -X-2/2/2 -Y-2/2/2 -Ostsyn.nc

# 剖分尺寸须同时为零或同时为正数
expect_fail grt static syn -Gstgrn_multi.nc -Cfaults.inp+i1/0 -Qrcv_points.txt -Ostsyn.nc

# 未指定点源源强
expect_fail grt static syn -Gstgrn.nc -Ostsyn.nc
# 使用接收点文件时，不可再设置接收深度
expect_fail grt static syn -Gstgrn_multi.nc -S1e20 -Ds2 -Dr0 -Qrcv_points.txt -Ostsyn.nc
# 接收断层的剖分尺寸须同时为零或同时为正数
expect_fail grt static syn -Gstgrn_multi.nc -S1e20 -Ds2 -Urcv_faults.inr+i0/1 -Ostsyn.nc

python -u test_static_syn.py

rm -rf stgrn.nc stgrn_multi.nc stsyn.nc
remove_test_files
