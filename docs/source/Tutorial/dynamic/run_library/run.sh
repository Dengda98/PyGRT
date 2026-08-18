#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

rm -rf GRN syn_root syn_subdir

# BEGIN GRN
# -Ds 和 -Dr 传入多个深度，程序会遍历全部深度组合
grt greenfn -Mmilrow -Ds2,4 -Dr0,2 -N256/0.02 -OGRN -R5,8,10
# END GRN

# BEGIN SYN ROOT
# 根目录模式下按震源深度、台站深度和震中距精确选择子目录
grt syn -GGRN -Ds4 -Dr2 -R8 -S1e24 -A30 -Osyn_root
# END SYN ROOT

# BEGIN SYN SUBDIR
# 子目录已经指定了三个几何参数，不能再设置 -Ds/-Dr/-R
grt syn -GGRN/milrow_4_2_8 -S1e24 -A30 -Osyn_subdir
# END SYN SUBDIR

rm -rf GRN syn_root syn_subdir
