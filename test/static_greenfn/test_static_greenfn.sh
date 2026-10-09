#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt static greenfn -h

# 积分控制、导数、统计、边界条件
grt static greenfn -M../milrow -D2/0 -R0/4/2 -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/2 -e -L10 -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/2 -K+k4+f+e1e-3 -S -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/2 -Cd -BrF -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/2 -Cp -BhR -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/2 -Cn -BrH -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R100,120 -L+l10+o2 -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R100,120 -L+a1e-3 -Ostgrn.nc

# 网格、距离、多个深度
grt static greenfn -M../milrow -D2/0 -X-2/2/2 -Y-2/2/2 -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R2,5,10 -Ostgrn.nc
grt static greenfn -M../milrow -Ds1,2 -Dr0/1/1 -Rdists -Ostgrn.nc
grt static greenfn -M../milrow -Dsdepsrc -Drdeprcv -R0/4/2 -e -Ostgrn.nc

# 不可同时设置单深度和多深度参数
expect_fail grt static greenfn -M../milrow -D2/0 -Ds1,2 -Dr0 -R0/4/2 -Ostgrn.nc
# 未指定多深度参数对应的接收深度
expect_fail grt static greenfn -M../milrow -Ds1,2 -R0/4/2 -Ostgrn.nc
# 震源深度不可为负数
expect_fail grt static greenfn -M../milrow -D-1/0 -R0/4/2 -Ostgrn.nc
# 震中距列表不可包含逆序或重复值
expect_fail grt static greenfn -M../milrow -D2/0 -R5,2 -Ostgrn.nc

# 接收网格不可只设置北向坐标
expect_fail grt static greenfn -M../milrow -D2/0 -X-2/2/2 -Ostgrn.nc

python -u test_static_greenfn.py

rm -rf stgrn.nc stgrtstats
remove_test_files
