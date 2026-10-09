#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt travt -h
grt travt -M../milrow -D2/0 -R5
grt travt -M../milrow -D2/0 -R2,5,10
grt travt -M../milrow -D2/0 -R2/10/4
grt travt -M../milrow -D2/0 -Rdists

# 震中距列表不可包含逆序或重复值
expect_fail grt travt -M../milrow -D2/0 -R5,2
# 震源深度不可为负数
expect_fail grt travt -M../milrow -D-1/0 -R5
# 未指定模型文件
expect_fail grt travt -D2/0 -R5

python -u test_travt.py

remove_test_files
