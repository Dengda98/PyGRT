#!/bin/bash

set -euo pipefail

rm -rf GRN* syn* *.svg *.tar.gz

# BEGIN
dist=30
depsrc=3

nt=1500
dt=0.05

modname="mod_sea"
out="GRN_sea"

# =================================================
# compute displacement in solid
deprcv=1.001
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist}

# =================================================
# compute displacement and its derivatives in liquid
deprcv=0.999
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -e
# just give a test focal mechanism
grt syn -G${out}/${modname}_${depsrc}_${deprcv}_${dist} -A40 -S1e20 -M100/55/120 -Osyn -e
# compute stress
grt stress syn

# END

python plot.py

cp stress.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf GRN* syn*
