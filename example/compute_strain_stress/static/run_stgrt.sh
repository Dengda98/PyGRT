#!/bin/bash

depsrc=2
deprcv=0

x1=-3
x2=3
nx=10

y1=-3
y2=3
ny=10

stgrt -Mmilrow -D${depsrc}/${deprcv} -X$x1/$x2/$nx -Y$y1/$y2/$ny -e > grn 

# Fault
S="u1e8"
stk=77
dip=88
rak=99
stgrt.syn -S$S -M$stk/$dip/$rak -e -N < grn > syn 

stgrt.strain < syn > strain
stgrt.stress < syn > stress