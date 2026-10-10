#!/bin/bash

set -euo pipefail
source ../common.sh

grt lamb2 -h
grt lamb2 -P0.25 -T0/2/0.05 -R10 -Ds5 -A30 -S+ssource.txt+rreceiver.txt+mmixed.txt > response.txt
grt lamb2 -P0.25 -T0/2/0.05 -R10 -Dr5 -A30 > response.txt
grt lamb2 -P0.25 -T0/2/0.05 -R10 -Ds5 -A30 -LP,S,SP > response.txt
grt lamb2 -P0.25 -T0/2/0.05 -R10 -Dr5 -A30 -LP,S,PS > response.txt

# 第二类 Lamb 问题不可同时设置地下源和地下接收点
expect_fail grt lamb2 -P0.25 -T0/2/0.05 -R10 -Ds5 -Dr3 -A30
# 地下源的深度不可为零或负数
expect_fail grt lamb2 -P0.25 -T0/2/0.05 -R10 -Ds0 -A30
# 水平距离不可为零或负数
expect_fail grt lamb2 -P0.25 -T0/2/0.05 -R-1 -Ds5 -A30

python -u test_lamb2.py

rm -rf response.txt source.txt receiver.txt mixed.txt
