#!/bin/bash
set -euo pipefail

# BEGIN GRN
# 不同震源深度、接收点深度和震中距的格林函数会在
# GRN/milrow_{depsrc}_{deprcv}_{dist}/ 路径下，使用SAC格式保存
grt greenfn -Mmilrow -D2/0 -N500/0.02 -OGRN -R5,8,10
# END GRN

# BEGIN grt.travt
grt travt -Mmilrow -D2/0 -R5,8,10
# END grt.travt

# BEGIN grt.sac2asc
grt sac2asc GRN/milrow_2_0_10/HFZ.sac > HFZ
# END grt.sac2asc

grt travt -Mmilrow -D2/0 -R5,8,10 > travt
head -n 10 HFZ > HFZ_head
echo "..." >> HFZ_head

python run.py

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -rf GRN HFZ
