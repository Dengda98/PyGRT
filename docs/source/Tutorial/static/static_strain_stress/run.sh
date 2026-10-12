#!/bin/bash

set -euo pipefail

# BEGIN
depsrc=2
deprcv=0

x1=-3
x2=3
dx=0.15

y1=-2.5
y2=2.5
dy=0.15

# -R 直接指定震中距，-e 表示计算空间导数
grt static greenfn -Mmilrow -D${depsrc}/${deprcv} -R0/10/0.15 -e -Ostgrn.nc

# -N 表示输出 ZNE 分量，-X/-Y 指定二维接收网格
grt static syn -S1e24 -M33/50/120 -e -N -Gstgrn.nc -X$x1/$x2/$dx -Y$y1/$y2/$dy -Ostsyn_dc_zne.nc

# 计算应变
grt static strain -Gstsyn_dc_zne.nc

# 计算旋转
grt static rotation -Gstsyn_dc_zne.nc

# 计算应力
grt static stress -Gstsyn_dc_zne.nc
# END

# BEGIN SPROJ
# 接收断层形态与震源机制独立设置
grt static sproj -Gstsyn_dc_zne.nc -M90/60/0
# END SPROJ

# BEGIN COULOMB
grt static coulomb -Gstsyn_dc_zne.nc -F0.4
# END COULOMB

python run.py

python plot.py

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -f stgrn.nc stsyn_*.nc
