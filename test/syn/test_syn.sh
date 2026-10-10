#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

dt=0.1

grt syn -h
grt greenfn -M../milrow -D2/0 -N32/"$dt" -R5 -e -OGRN -s

# 点源机制、源强、时间函数、积分和微分
grt syn -GGRN -A22 -S1e20 -Osyn
grt syn -GGRN/milrow_2_0_5 -A22 -S1e16 -F-1/2/-4 -Osyn
grt syn -GGRN -A22 -S1e20 -M33/44/55 -Osyn
grt syn -GGRN -A22 -Su1e10 -M33/44 -Osyn
grt syn -GGRN -A22 -S1e20 -T1/-2/-5/0.5/3/1.2 -Osyn
grt syn -GGRN -A22 -S1e20 -Dp/0.6 -Osyn
grt syn -GGRN -A22 -S1e20 -Dt/0.2/0.2/0.3 -Osyn
grt syn -GGRN -A22 -S1e20 -Dc/0.2/0.4 -Osyn
grt syn -GGRN -A22 -S1e20 -Dr/1.2 -Osyn
# 时间函数按格林函数的 dt 采样
cat > time_function.txt <<'EOF'
0.0
2.5
5.0
2.5
0.0
EOF
grt syn -GGRN -A22 -S1e20 -D0/time_function.txt -Osyn
# 两列输入从 0.0 开始，支持不同采样间隔和非等距采样
cat > time_function_times.txt <<'EOF'
0.00  0.0
0.05  1.0
0.10  2.0
0.15  3.0
0.20  4.0
0.25  3.0
0.30  2.0
0.35  1.0
0.40  0.0
EOF
cat > time_function_nonuniform.txt <<'EOF'
# time amplitude

0.00  0.0
0.07  1.0
0.23  3.0
0.45  0.0
EOF
grt syn -GGRN -A22 -S1e20 -D0/time_function_times.txt -Osyn
grt syn -GGRN -A22 -S1e20 -D0/time_function_nonuniform.txt+d0.2 -Osyn

grt syn -GGRN -A22 -S1e20 -I1 -Osyn
grt syn -GGRN -A22 -S1e20 -J1 -N -e -Osyn

# 多深度库查询、插值、任意接收点和有限断层
grt greenfn -M../milrow -Ds1,2 -Dr0,1 -N16/"$dt" -R1,6,12 -e -OGRN_MULTI -s
grt syn -GGRN_MULTI -Ds1.5 -Dr0.5 -R8 -A22 -S1e20 -Osyn
grt syn -GGRN_MULTI -Ds1.5 -Dr0.5 -R8 -A22 -S1e20 -i0 -Osyn
grt syn -GGRN_MULTI -Ds2 -S1e20 -Qrcv_points.txt -P2 -Osyn_points
grt syn -GGRN_MULTI -Ds2 -S1e20 -Qrcv_geometry.txt -e -Osyn_points
grt syn -GGRN_MULTI -Ds2 -S1e20 -Urcv_faults.inr -Osyn_faults
grt syn -GGRN_MULTI -Cfaults.inp+i1/1 -Qrcv_points.txt -e -Osyn_points
grt syn -GGRN_MULTI -Cfaults.inr -Urcv_faults.inr+i1/1 -Osyn_faults

# 各断层的破裂过程和延迟，以及全局时间函数覆盖行末设置
grt syn -GGRN_MULTI -Cfaults_rupture.inr+i1/1 -Qrcv_points.txt -Osyn_points
grt syn -GGRN_MULTI -Cfaults_rupture.inr+i1/1 -Qrcv_points.txt -Dp/0.4+d0.1 -Osyn_points
grt syn -GGRN_MULTI -Cfaults_rupture.inr+i1/1 -Qrcv_points.txt -D0/time_function_nonuniform.txt -Osyn_points

# 有限断层的时间函数不可只设置在部分记录中
expect_fail grt syn -GGRN_MULTI -Cfaults_partial.inr+i1/1 -Qrcv_points.txt -Osyn_points
# 破裂延迟不可为负数
expect_fail grt syn -GGRN_MULTI -Cfaults_negative_delay.inr+i1/1 -Qrcv_points.txt -Osyn_points

# 未指定多深度、多距离库对应的震源深度、接收深度和距离
expect_fail grt syn -GGRN_MULTI -A22 -S1e20 -Osyn
# 震中距不可超出格林函数库范围
expect_fail grt syn -GGRN_MULTI -Ds2 -Dr0 -R13 -A22 -S1e20 -Osyn
# 不可同时设置单力源和双力偶源机制
expect_fail grt syn -GGRN -A22 -S1e20 -F1/2/3 -M33/44/55 -Osyn
# 不可同时设置接收点文件和接收断层文件
expect_fail grt syn -GGRN_MULTI -Ds2 -S1e20 -Qrcv_points.txt -Urcv_faults.inr -Osyn

# 未指定点源源强
expect_fail grt syn -GGRN -A22 -Osyn
# 使用接收点文件时，不可再设置接收深度
expect_fail grt syn -GGRN_MULTI -Ds2 -Dr0 -S1e20 -Qrcv_points.txt -Osyn
# 使用格林函数节点子目录时，不可再设置深度和距离选择参数
expect_fail grt syn -GGRN/milrow_2_0_5 -Ds2 -A22 -S1e20 -Osyn
# 抛物线时间函数的持续时间不可为零或负数
expect_fail grt syn -GGRN -A22 -S1e20 -Dp/0 -Osyn

python -u test_syn.py

rm -rf GRN GRN_MULTI syn syn_points syn_faults time_function.txt
rm -f time_function_times.txt time_function_nonuniform.txt
remove_test_files
