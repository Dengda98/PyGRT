#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt strain -h
grt stress -h
grt rotation -h
grt sproj -h
grt coulomb -h

grt greenfn -M../milrow -Ds2 -Dr0,1 -N16/0.1 -R1,6,12 -e -OGRN -s

grt syn -GGRN -Ds2 -Dr0 -R6 -A22 -S1e20 -Osyn_plain
# 未计算位移导数，不可计算应变
expect_fail grt strain -Gsyn_plain

# ZRT 和 ZNE 单台后处理
grt syn -GGRN -Ds2 -Dr0 -R6 -A22 -S1e20 -e -Osyn
grt strain -Gsyn
grt stress -Gsyn
grt rotation -Gsyn
# 未指定接收断层几何
expect_fail grt sproj -Gsyn
# 未计算法向和剪切应力投影，不可计算 Coulomb 应力
expect_fail grt coulomb -Gsyn -F0.4
grt sproj -Gsyn -M33/44/55
grt coulomb -Gsyn -F0.4

grt syn -GGRN -Ds2 -Dr0 -R6 -A22 -S1e20 -N -e -Osyn
grt strain -Gsyn
grt stress -Gsyn
grt rotation -Gsyn
grt sproj -Gsyn -M33/44/55
grt coulomb -Gsyn -F0

# 多台站文件中的几何和外部几何
grt syn -GGRN -Ds2 -S1e20 -Qrcv_geometry.txt -e -Osyn_points
grt strain -Gsyn_points
grt stress -Gsyn_points
grt rotation -Gsyn_points
grt sproj -Gsyn_points
grt sproj -Gsyn_points -Qrcv_geometry.txt
grt coulomb -Gsyn_points -F0.4

# 有限接收断层的滑动角覆盖
grt syn -GGRN -Ds2 -S1e20 -Urcv_faults.inr -e -Osyn_faults
grt stress -Gsyn_faults
grt sproj -Gsyn_faults
grt sproj -Gsyn_faults -M55+f
grt coulomb -Gsyn_faults -F0.4

# 未指定动态合成目录
expect_fail grt strain
# 输入的动态合成目录不存在
expect_fail grt stress -Gmissing
# 不可设置不支持的选项 -x
expect_fail grt rotation -Gsyn -x
# 未指定接收断层的滑动角
expect_fail grt sproj -Gsyn -M33/44
# 摩擦系数不可为负数
expect_fail grt coulomb -Gsyn -F-0.1

rm -rf syn syn_points syn_faults

python -u test_tensors.py

rm -rf GRN syn syn_points syn_faults syn_plain
remove_test_files
