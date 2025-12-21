#!/bin/bash

set -euo pipefail

rm -rf GRN* syn* *.svg *.tar.gz


dist=10
deprcv=0.5
depsrc=3.5

nt=1500
dt=0.01

modname="water"
out="GRN"

grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -Ge
grt syn -G${out}/${modname}_${depsrc}_${deprcv}_${dist} -S1e20 -A0 -I1 -Osyn

python plot.py

cp disp.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf GRN* syn*
