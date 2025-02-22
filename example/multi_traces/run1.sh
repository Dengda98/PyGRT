#!/bin/bash 

set -eu

# 震中距数组
distarr=($(seq -f "%.2f" 0.5 0.01 1.5  | tr '\n' ','))
distarr=${distarr%,}  # 删除最后的逗号

# 震源台站都在地表
depsrc=0
deprcv=0

# 时窗
nt=3000
dt=0.02  # 50Hz

modname="mod1"
out="GRN"

# compute Green's Functions
grt -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${distarr} -G0/1/0/0  # 只生成垂直力源的格林函数
