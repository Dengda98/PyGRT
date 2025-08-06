#!/bin/bash
set -eu

rm -rf grn* syn* strain stress rotation *.png

depsrc=2
deprcv=0

x1=-3
x2=3
dx=0.15

y1=-2.5
y2=2.5
dy=0.15
# -e 表示计算空间导数
grt static greenfn -Mmilrow -D${depsrc}/${deprcv} -X$x1/$x2/$dx -Y$y1/$y2/$dy -e > grn

# -N 表示输出ZNE分量
grt static syn -S1e24 -M33/50/120 -e -N < grn > syn_dc_zne

# 计算应变
grt static strain < syn_dc_zne > strain

# 计算旋转
grt static rotation < syn_dc_zne > rotation

# 计算应力
grt static stress < syn_dc_zne > stress