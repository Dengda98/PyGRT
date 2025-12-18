#!/bin/bash

rm -rf GRN*

modname="mod"

cat > ${modname} <<EOF
0   8   4.62  3.3
EOF

dist=10
depsrc=0
deprcv=0

nt=1000
dt=0.005

out="GRN"

rm GRN* -rf

zeta=0.8

L="-L20"
# 取不同震源深度，输出核函数的变化 seq 0 0.02 0.2; 
for depsrc in $(seq 0 0.02 0.4); do
    grt greenfn -M${modname} -O${out} -N${nt}/${dt}/$zeta -D${depsrc}/${deprcv} -R${dist} $L -H5/5 -S25 -K+k50+s5 -Cn
done

python plot_depth_kernel.py