#!/bin/bash 

depsrc=2
nt=1024
dt=0.01

modname="milrow"
out="GRN"


diff=0.001
deprcv=3.2
deprcv2=`echo $deprcv+$diff | bc | awk '{printf("%.3f", $1)}'`
dist=10
dist2=`echo $dist+$diff | bc | awk '{printf("%.3f", $1)}'`

grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/$deprcv -R$dist -e
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/$deprcv2 -R$dist
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/$deprcv -R$dist2
