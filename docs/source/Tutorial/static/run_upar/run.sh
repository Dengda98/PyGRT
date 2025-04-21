#!/bin/bash
set -eu

depsrc=2
deprcv=0

x1=-3
x2=3
nx=41

y1=-2.5
y2=2.5
ny=33
# -e 表示计算空间导数
stgrt -Mmilrow -D${depsrc}/${deprcv} -X$x1/$x2/$nx -Y$y1/$y2/$ny -e > grn

# -N 表示输出ZNE分量
stgrt.syn -S1e24 -M33/50/120 -e -N < grn > syn_dc_zne

# 计算应变
stgrt.strain < syn_dc_zne > strain

# 计算旋转
stgrt.rotation < syn_dc_zne > rotation

# 计算应力
stgrt.stress < syn_dc_zne > stress