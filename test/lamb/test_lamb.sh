#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

dt=0.1

grt lamb -h

# 三类源点/接收点深度组合，单力、双力偶、矩张量和张裂源
grt lamb -H8/4.62/3.3 -N32/"$dt" -R10 -Ds0 -Dr0 -A30 -F0.5/-1/2 -S1e15 -Olamb
grt lamb -H8/4.62/3.3 -N32/"$dt" -R10 -Ds5 -Dr0 -A30 -M100/30/70 -S1e20 -e -n -Olamb
grt lamb -H8/4.62/3.3 -N32/"$dt" -R10 -Ds0 -Dr5 -A30 -T1/0/0/1/0/1 -S1e20 -Olamb
grt lamb -H8/4.62/3.3 -N32/"$dt" -R10 -Ds5 -Dr1 -A30 -M100/30 -Su1e10 -Olamb

# 相位、时间函数、积分微分、时间偏移
grt lamb -H8/4.62/3.3 -N32/"$dt" -R10 -Ds5 -Dr0 -A30 -S1e20 -LP,S,SP -Olamb
grt lamb -H8/4.62/3.3 -N32/"$dt" -R10 -Ds5 -Dr1 -A30 -S1e20 -Dp/0.5 -I1 -J1 -E0.2/9 -Olamb
grt lamb -H8/4.62/3.3 -N32/"$dt" -R10 -Ds5 -Dr1 -A30 -S1e20 -Ep-1 -Olamb
# 自定义时间函数按 dt 采样
cat > time_function.txt <<'EOF'
0.0
2.5
5.0
2.5
0.0
EOF
grt lamb -H8/4.62/3.3 -N32/"$dt" -R10 -Ds5 -Dr1 -A30 -S1e20 -D0/time_function.txt+d0.2 -Olamb
# 两列时间、振幅输入支持非等距采样
cat > time_function_nonuniform.txt <<'EOF'
0.00  0.0
0.07  1.0
0.23  3.0
0.45  0.0
EOF
grt lamb -H8/4.62/3.3 -N32/"$dt" -R10 -Ds5 -Dr1 -A30 -S1e20 -D0/time_function_nonuniform.txt+d0.2 -Olamb

# 任意接收点、有限震源和接收断层
grt lamb -H8/4.62/3.3 -N16/"$dt" -Ds5 -S1e20 -Qrcv_geometry.txt -e -P2 -Olamb_points
grt lamb -H8/4.62/3.3 -N16/"$dt" -Ds5 -S1e20 -Urcv_faults.inr+i1/1 -Olamb_faults
grt lamb -H8/4.62/3.3 -N16/"$dt" -Cfaults.inp+i1/1 -Qrcv_points.txt -Olamb_points -s
grt lamb -H8/4.62/3.3 -N16/"$dt" -Cfaults.inr+i1/1 -Urcv_faults.inr -Olamb_faults -s

# 显式不再剖分，每条矩形源及接收断层只取一个中心点
grt lamb -H8/4.62/3.3 -N16/"$dt" -Cfaults.inp+i0/0 -Urcv_faults.inr+i0/0 -Olamb_faults -s

# 各断层的破裂过程和延迟，以及全局时间函数覆盖行末设置
grt lamb -H8/4.62/3.3 -N16/"$dt" -Cfaults_rupture.inr+i1/1 -Qrcv_points.txt -Olamb_points -s
grt lamb -H8/4.62/3.3 -N16/"$dt" -Cfaults_rupture.inr+i1/1 -Qrcv_points.txt -Dt/0.2/0.1/0.3+d0.1 -Olamb_points -s

# 有限断层的时间函数不可只设置在部分记录中
expect_fail grt lamb -H8/4.62/3.3 -N16/"$dt" -Cfaults_partial.inr+i1/1 -Qrcv_points.txt -Olamb_points
# 破裂延迟不可为负数
expect_fail grt lamb -H8/4.62/3.3 -N16/"$dt" -Cfaults_negative_delay.inr+i1/1 -Qrcv_points.txt -Olamb_points

# 源和接收点都在地表时，不可使用爆炸源
expect_fail grt lamb -H8/4.62/3.3 -N16/"$dt" -R10 -Ds0 -Dr0 -A30 -S1e20 -Olamb
# 水平距离不可为零或负数
expect_fail grt lamb -H8/4.62/3.3 -N16/"$dt" -R0 -Ds5 -Dr1 -A30 -S1e20 -Olamb
# 不可同时设置接收点文件和接收断层文件
expect_fail grt lamb -H8/4.62/3.3 -N16/"$dt" -Ds5 -S1e20 -Qrcv_points.txt -Urcv_faults.inr -Olamb

# 剖分尺寸不可为负数
expect_fail grt lamb -H8/4.62/3.3 -N16/"$dt" -Cfaults.inp+i-1/-1 -Qrcv_points.txt -Olamb_points

# 未指定点源源强
expect_fail grt lamb -H8/4.62/3.3 -N16/"$dt" -R10 -Ds5 -Dr1 -A30 -Olamb
# 矩形有限震源不可省略剖分尺寸
expect_fail grt lamb -H8/4.62/3.3 -N16/"$dt" -Cfaults.inr -Qrcv_points.txt -Olamb_points

python -u test_lamb.py

rm -rf lamb lamb_points lamb_faults time_function.txt time_function_nonuniform.txt
remove_test_files
