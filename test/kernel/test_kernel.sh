#!/bin/bash

set -euo pipefail
source ../common.sh

grt kernel -h
grt ker2asc -h
grt kernel -M../thin1 -D0.03/0 -F0/1/0.5 -C0.1 -OKERNEL
grt kernel -M../thin1 -D0.03/0 -F0/1/0.5+w1 -C0.2/2/0.2 -BrH -e -P2 -OKERNEL
grt kernel -M../seafloor -D5/0.1 -F0/0.1/0.05 -C0.5 -OKERNEL

# ker2asc 读取波数积分统计文件
grt greenfn -M../milrow -D2/0 -N8/0.1 -R5 -S1 -OGRN -s
for file in GRN_grtstats/milrow_2_0/*; do
    grt ker2asc "$file" > kernel.txt
done

# 频率采样间隔不可为零或负数
expect_fail grt kernel -M../thin1 -D0.03/0 -F0/1/0 -C0.1 -OKERNEL
# 相速度采样间隔不可为零或负数
expect_fail grt kernel -M../thin1 -D0.03/0 -F0/1/0.5 -C0 -OKERNEL
# 输入的积分统计文件不存在
expect_fail grt ker2asc missing

rm -rf KERNEL GRN GRN_grtstats kernel.txt
