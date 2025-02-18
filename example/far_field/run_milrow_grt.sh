#!/bin/bash


dist=1800
depsrc=2
deprcv=0

nt=1400
dt=1

T0=190

modname="milrow"
out="GRN"

grt -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -H0.004/-1 -L20 -E$T0

out="GRN_filon"  
grt -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -H0.004/-1 -L-5 -E$T0 # Filon's Integration
