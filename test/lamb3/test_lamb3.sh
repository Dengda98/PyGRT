#!/bin/bash

set -euo pipefail
source ../common.sh

grt lamb3 -h
grt lamb3 -P0.25 -T0/2/0.05 -R10 -D2/1 -A30 -S+ssource.txt+rreceiver.txt+mmixed.txt > response.txt
grt lamb3 -P0.25 -T0/2/0.05 -R10 -D2/1 -A30 -LP,S,PP,SS,PS,SP,sPs > response.txt

# 地下源的深度不可为零或负数
expect_fail grt lamb3 -P0.25 -T0/2/0.05 -R10 -D0/1 -A30
# 地下接收点的深度不可为零或负数
expect_fail grt lamb3 -P0.25 -T0/2/0.05 -R10 -D2/0 -A30
# 水平距离不可为零或负数
expect_fail grt lamb3 -P0.25 -T0/2/0.05 -R0 -D2/1 -A30

python -u test_lamb3.py

rm -rf response.txt source.txt receiver.txt mixed.txt
