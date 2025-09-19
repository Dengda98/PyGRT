#!/bin/bash

depsrc=2
deprcv=0

x1=-3
x2=3
dx=0.6

y1=-3
y2=3
dy=0.6

grt static greenfn -Mmilrow -D${depsrc}/${deprcv} -X$x1/$x2/$dx -Y$y1/$y2/$dy -e -Ostgrn.nc 

# Fault
S="u1e8"
stk=77
dip=88
rak=99
grt static syn -S$S -M$stk/$dip/$rak -e -N -Gstgrn.nc -Ostsyn.nc 

grt static strain stsyn.nc 
grt static stress stsyn.nc 