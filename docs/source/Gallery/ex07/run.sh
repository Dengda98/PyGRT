#!/bin/bash 

set -euo pipefail

rm -rf GRN* *.svg *.tar.gz

# 震中距数组
# 输出到文件
seq -f "%.2f" 0.5 0.01 1.5 > dists
distarr=($(cat dists | tr '\n' ','))
distarr=${distarr%,}  # 删除最后的逗号

# 震源台站都在地表
depsrc=0
deprcv=0

# 时窗
nt=800
dt=0.02  # 50Hz

modname="mod"
out="GRN"

# compute Green's Functions
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -Rdists -Gv  # 只生成垂直力源的格林函数

python plot.py

cp multi_traces.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf GRN* dists