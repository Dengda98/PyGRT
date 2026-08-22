#!/bin/bash

set -euo pipefail

rm -rf stgrn* stsyn* *.svg

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
grt static strain stsyn_dc_zne.nc

# 计算旋转
grt static rotation stsyn_dc_zne.nc

# 计算应力
grt static stress stsyn_dc_zne.nc
# END

# C 示例与 Python 示例共用文件名，跑 Python 前清掉残留
rm -rf stgrn* stsyn*
python run.py

rm -rf stgrn* stsyn*
