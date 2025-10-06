#!/bin/bash 

# 震中距数组
# 输出到文件
seq -f "%.2f" 0.5 0.01 1.5 > dists
distarr=($(cat dists | tr '\n' ','))
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
# OPTION-1: input distances array
# grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${distarr} -Gv  # 只生成垂直力源的格林函数
# OPTION-2: 
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -Rdists -Gv  # 只生成垂直力源的格林函数
