#!/bin/bash

set -euo pipefail
source ../common.sh

grt lamb1 -h
grt lamb1 -P0.25 -T0/2/0.05 -A30 > response.txt
grt lamb1 -P0.25 -T0/2/0.05 -A30 -LP,S > response.txt
grt lamb1 -P0.25 -T0/2/0.05 -A30 -LR > response.txt
grt lamb1 -P0.25 -T0/2/0.05 -A30 -C2e-4 > response.txt

# 移动源的无量纲速度不可为零或负数
expect_fail grt lamb1 -P0.25 -T0/2/0.05 -A30 -C0
# 无量纲时间采样间隔不可为零或负数
expect_fail grt lamb1 -P0.25 -T0/2/0 -A30

# 泊松比不可超出 (0, 0.5) 范围
expect_fail grt lamb1 -P0.5 -T0/2/0.05 -A30
# 移动源的接收点不可位于 x1 轴上
expect_fail grt lamb1 -P0.25 -T0/2/0.05 -A0 -C2e-4

python -u test_lamb1.py

rm -rf response.txt
