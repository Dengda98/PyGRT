#!/bin/bash

set -euo pipefail

# -----------------------------------------------------------------------------------
# BEGIN CLI RCVFN
# P 波使用水平射线参数，SV 波使用底部半空间的入射角
grt rcvfn -Mmod1 -P0.03 -TP -N500/0.1 -A5.0 -E-10 -W -OCLI_P
grt rcvfn -Mmod1 -I10/0 -TS -N500/0.1 -A5.0 -E-10 -W -OCLI_S
# END CLI RCVFN
# -----------------------------------------------------------------------------------

python run.py

python plot.py

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -rf CLI_* PY_*
