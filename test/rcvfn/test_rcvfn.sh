#!/bin/bash

set -euo pipefail
source ../common.sh

grt rcvfn -h
grt rcvfn -M../milrow -P0.12 -TP -N32/0.1 -Oreceiver
grt rcvfn -M../milrow -I20/2 -TS -N32/0.1+w0.9+n2+a+f -A0.8 -E0.5 -W -Oreceiver

# 不可同时设置水平射线参数和入射角
expect_fail grt rcvfn -M../milrow -P0.12 -I20 -TP -N32/0.1 -Oreceiver
# 未指定水平射线参数或入射角
expect_fail grt rcvfn -M../milrow -TP -N32/0.1 -Oreceiver
# 入射角不可达到或超过 90 度
expect_fail grt rcvfn -M../milrow -I90 -TP -N32/0.1 -Oreceiver

# 水平射线参数不可大到使 P 波无法传播
expect_fail grt rcvfn -M../milrow -P1 -TP -N32/0.1 -Oreceiver

python -u test_rcvfn.py

rm -rf receiver
