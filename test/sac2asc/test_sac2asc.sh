#!/bin/bash

set -euo pipefail
source ../common.sh

grt sac2asc -h
grt greenfn -M../milrow -D2/0 -N16/0.1 -R5 -OGRN -s
grt sac2asc GRN/milrow_2_0_5/VFZ.sac > wave.txt
# 输入的 SAC 文件不存在
expect_fail grt sac2asc missing.sac

rm -rf GRN wave.txt
