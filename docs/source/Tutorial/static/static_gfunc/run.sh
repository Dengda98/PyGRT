#!/bin/bash
set -euo pipefail

# BEGIN GRN
grt static greenfn -Mmilrow -D2/0 -R0/10/0.1 -Ostgrn.nc
# END GRN

ncdump -h stgrn.nc > grn_head

python run.py

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -f stgrn.nc
