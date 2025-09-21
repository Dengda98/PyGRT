#!/bin/bash

dist=10
depsrc=3

nt=2000
dt=0.01

modname="mod_sea"
out="GRN_sea"

# compute displacement in solid
deprcv=1.001
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist}

# compute displacement and its derivatives in liquid
deprcv=0.999
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -e
# just give a random focal mechanism
grt syn -G${out}/${modname}_${depsrc}_${deprcv}_${dist} -A20 -S1e20 -M100/20/120 -OSYN -e
# compute stress
grt stress SYN
