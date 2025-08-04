#!/bin/bash

depsrc=2
deprcv=0

x1=-3
x2=3
dx=0.6

y1=-3
y2=3
dy=0.6

grt static greenfn -Mmilrow -D${depsrc}/${deprcv} -X$x1/$x2/$dx -Y$y1/$y2/$dy -e > grn 

# Fault
S="u1e8"
stk=77
dip=88
rak=99
grt static syn -S$S -M$stk/$dip/$rak -e -N < grn > syn 

grt static strain < syn > strain
grt static stress < syn > stress