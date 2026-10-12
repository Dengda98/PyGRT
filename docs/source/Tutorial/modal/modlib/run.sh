#!/bin/bash
set -euo pipefail

# BEGIN LIBRARY
# eigenv 的 -F 使用等间隔频率，供 modsum 进行逆傅里叶变换
grt eigenv -Mmilrow -F0/1/0.02 -SR -N -Cphase_R.nc
grt eigenv -Mmilrow -F0/1/0.02 -SL -N -Cphase_L.nc

# 计算全部震源深度、台站深度和震中距的组合
grt modsum -Cphase_R.nc -Ds2,4 -Dr0,2 -R80,100,120 -N0 -OGRN -W4
grt modsum -Cphase_L.nc -Ds2,4 -Dr0,2 -R80,100,120 -N0 -OGRN -W4
# END LIBRARY

python run.py

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -rf GRN* phase_*.nc
