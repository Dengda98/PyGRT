#!/bin/bash

set -euo pipefail

rm -rf GRN* syn* *.svg

dist=10
depsrc=2
deprcv=0

nt=1000
dt=0.01

modname="milrow"
out="GRN"

grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -e 

# Fault
S="u1e10"
stk=77
dip=88
rak=99
az=30
grt syn -G${out}/${modname}_${depsrc}_${deprcv}_${dist} -A$az -S$S -M$stk/$dip/$rak -Dp/0.6 -Osyn_dc -e

grt strain syn_dc
grt stress syn_dc
grt rotation syn_dc

python plot.py

cp stress.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf GRN* syn*