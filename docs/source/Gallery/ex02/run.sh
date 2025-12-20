#!/bin/bash 

set -euo pipefail

rm -rf *.svg *.tar.gz GRN*

depsrc=2
nt=1500
dt=0.01

modname="milrow"
out="GRN"


diff=0.001
deprcv=3.2
deprcv2=`echo $deprcv+$diff | bc | awk '{printf("%.3f", $1)}'`
dist=10
dist2=`echo $dist+$diff | bc | awk '{printf("%.3f", $1)}'`

grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/$deprcv -R$dist -e
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/$deprcv2 -R$dist
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/$deprcv -R$dist2

python plot.py

cp compare_z.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf GRN*