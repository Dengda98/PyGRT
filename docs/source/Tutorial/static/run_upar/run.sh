#!/bin/bash

# BEGIN
set -euo pipefail

rm -rf stgrn* stsyn* *.svg

depsrc=2
deprcv=0

x1=-3
x2=3
dx=0.15

y1=-2.5
y2=2.5
dy=0.15
# -e 表示计算空间导数
grt static greenfn -Mmilrow -D${depsrc}/${deprcv} -X$x1/$x2/$dx -Y$y1/$y2/$dy -e -Ostgrn.nc

# -N 表示输出ZNE分量
grt static syn -S1e24 -M33/50/120 -e -N -Gstgrn.nc -Ostsyn_dc_zne.nc

# 计算应变
grt static strain stsyn_dc_zne.nc

# 计算旋转
grt static rotation stsyn_dc_zne.nc

# 计算应力
grt static stress stsyn_dc_zne.nc
# END

python run.py

rm -rf stgrn* stsyn*