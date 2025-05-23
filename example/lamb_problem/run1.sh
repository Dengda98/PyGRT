#!/bin/bash 

dist=10

# both on surface
depsrc=0
deprcv=0

nt=500
dt=0.01

modname="halfspace"
out="GRN"

# compute Green's Functions
grt -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -G0/1/1/0   # Vertical and Horizontal Force
