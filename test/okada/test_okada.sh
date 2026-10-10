#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt okada -h

# 旧选项 -I 仍兼容，并提示改用 -H
grt okada -I6/3.464/2.7 -Su1e12 -Ds10 -Dr0 -X0/0/1 -Y0/0/1 -Ookada.nc -s

grt okada -H6/3.464/2.7 -Su1e12 -Ds10 -Dr0 -X-2/2/2 -Y-2/2/2 -Ookada.nc
grt okada -H6/3.464/2.7 -Su1e16 -Ds10 -Dr0 -X-2/2/2 -Y-2/2/2 -M100/20/80 -N -e -Ookada.nc
grt okada -H6/3.464/2.7 -Su1e16 -Ds10 -Dr0 -X-2/2/2 -Y-2/2/2 -M100/20 -Ookada.nc
grt okada -H6/3.464/2.7 -Su1e12 -Ds10 -Qrcv_geometry.txt -e -Ookada.nc
grt okada -H6/3.464/2.7 -Su1e12 -Ds10 -Urcv_faults.inr -Ookada.nc
grt okada -H6/3.464/2.7 -Su1e12 -Ds10 -Urcv_faults.inr+i0/0 -Ookada.nc
grt okada -H6/3.464/2.7 -Su1e12 -Ds10 -Urcv_faults.inr+i1/1 -Ookada.nc
grt okada -H6/3.464/2.7 -Cfaults.inp -Dr0 -X-2/2/2 -Y-2/2/2 -N -e -Ookada.nc
grt okada -H6/3.464/2.7 -Cfaults.inr -Qrcv_points.txt -Ookada.nc

# Okada 震源不可附加 +i，即使尺寸为零
expect_fail grt okada -H6/3.464/2.7 -Cfaults.inp+i0/0 -Qrcv_points.txt -Ookada.nc
expect_fail grt okada -H6/3.464/2.7 -Cfaults.inp+i1/1 -Qrcv_points.txt -Ookada.nc

# 震源深度不可为负数
expect_fail grt okada -H6/3.464/2.7 -Su1e12 -Ds-1 -Qrcv_points.txt -Ookada.nc
# 不可同时设置有限震源文件和点源源强
expect_fail grt okada -H6/3.464/2.7 -Su1e12 -Cfaults.inp -Qrcv_points.txt -Ookada.nc
# 不可同时设置接收点文件和接收断层文件
expect_fail grt okada -H6/3.464/2.7 -Su1e12 -Ds10 -Qrcv_points.txt -Urcv_faults.inr -Ookada.nc

# 接收网格不可只设置北向坐标
expect_fail grt okada -H6/3.464/2.7 -S1e20 -Ds10 -Dr0 -X-2/2/2 -Ookada.nc
# 使用接收点文件时，不可再设置接收深度
expect_fail grt okada -H6/3.464/2.7 -S1e20 -Ds10 -Dr0 -Qrcv_points.txt -Ookada.nc

python -u test_okada.py

rm -rf okada.nc
remove_test_files
