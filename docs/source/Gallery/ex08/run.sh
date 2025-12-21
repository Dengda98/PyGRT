#!/bin/bash 

set -euo pipefail

rm -rf GRN* syn* *.svg *.tar.gz

dist=3

depsrc=3
deprcv=0

nt=1000
dt=0.005

# ----------------------- model with soil ---------------------
modname="killari"
out="GRN"

# compute Green's Functions
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist}

# generate synthetic seismograms, integrate once
grt syn -G${out}/${modname}_${depsrc}_${deprcv}_${dist} -A45 -S1e24 -M0/45/90 -I1 -Osyn


# ----------------------- model without soil ---------------------
modname="killari_nosoil"
out="GRN_nosoil"

# compute Green's Functions
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist}

# generate synthetic seismograms, integrate once
grt syn -G${out}/${modname}_${depsrc}_${deprcv}_${dist} -A45 -S1e24 -M0/45/90 -I1 -Osyn_nosoil

python plot.py

cp compare.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf GRN* syn*