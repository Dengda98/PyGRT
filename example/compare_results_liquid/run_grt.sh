#!/bin/bash

dist=100
depsrc=20
deprcv=8.2

nt=512
dt=0.2

modname="liquid"
out="GRN"

grt -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -L20

