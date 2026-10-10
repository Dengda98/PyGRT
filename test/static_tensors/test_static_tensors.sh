#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt static strain -h
grt static stress -h
grt static rotation -h
grt static sproj -h
grt static coulomb -h

grt static greenfn -M../milrow -Ds2 -Dr0,1 -R0/12/3 -e -Ostgrn.nc

grt static syn -Gstgrn.nc -S1e20 -Dr0 -Ostsyn_plain.nc
# 未计算位移导数，不可计算应变
expect_fail grt static strain stsyn_plain.nc

# 网格结果的 ZRT 和 ZNE 后处理
grt static syn -Gstgrn.nc -S1e20 -Dr0 -X-2/2/2 -Y-2/2/2 -e -Ostsyn.nc
grt static strain stsyn.nc
grt static stress stsyn.nc
grt static rotation stsyn.nc
# 未指定接收断层几何
expect_fail grt static sproj -Gstsyn.nc
# 未计算法向和剪切应力投影，不可计算 Coulomb 应力
expect_fail grt static coulomb -Gstsyn.nc -F0.4
grt static sproj -Gstsyn.nc -M33/44/55
grt static coulomb -Gstsyn.nc -F0.4

grt static syn -Gstgrn.nc -S1e20 -Dr0 -X-2/2/2 -Y-2/2/2 -N -e -Ostsyn.nc
grt static strain stsyn.nc
grt static stress stsyn.nc
grt static rotation stsyn.nc
grt static sproj -Gstsyn.nc -M33/44/55
grt static coulomb -Gstsyn.nc -F0

# 任意接收点中的几何和外部几何
grt static syn -Gstgrn.nc -S1e20 -Qrcv_geometry.txt -e -Ostsyn_points.nc
grt static strain stsyn_points.nc
grt static stress stsyn_points.nc
grt static rotation stsyn_points.nc
grt static sproj -Gstsyn_points.nc
grt static sproj -Gstsyn_points.nc -Qrcv_geometry.txt
grt static coulomb -Gstsyn_points.nc -F0.4

# 有限接收断层的滑动角覆盖
grt static syn -Gstgrn.nc -S1e20 -Urcv_faults.inr -e -Ostsyn_faults.nc
grt static stress stsyn_faults.nc
grt static sproj -Gstsyn_faults.nc
grt static sproj -Gstsyn_faults.nc -M55+f
grt static coulomb -Gstsyn_faults.nc -F0.4

# 输入的静态合成文件不存在
expect_fail grt static strain missing.nc
# 输入的静态合成文件不存在
expect_fail grt static stress missing.nc
# 输入的静态合成文件不存在
expect_fail grt static rotation missing.nc
# 未指定接收断层的滑动角
expect_fail grt static sproj -Gstsyn.nc -M33/44
# 摩擦系数不可为负数
expect_fail grt static coulomb -Gstsyn.nc -F-0.1

python -u test_static_tensors.py

rm -rf stgrn.nc stsyn.nc stsyn_points.nc stsyn_faults.nc stsyn_plain.nc
remove_test_files
