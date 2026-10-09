#!/bin/bash
set -euo pipefail

rm -rf GRN syn_linear syn_nearest

# BEGIN GRN
# 计算全部震源深度、台站深度和震中距的组合
grt greenfn -Mmilrow -Ds2,4 -Dr0,2 -N256/0.02 -OGRN -R5,8,10
# END GRN

# BEGIN SYN
# 目标位置位于库范围内，无需与采样节点重合
# -i1 使用线性插值，也是默认方式
grt syn -GGRN -Ds3 -Dr1 -R7 -S1e24 -A30 -i1 -Osyn_linear

# 同一目标位置，改用 -i0 选择最近邻
grt syn -GGRN -Ds3 -Dr1 -R7 -S1e24 -A30 -i0 -Osyn_nearest
# END SYN
