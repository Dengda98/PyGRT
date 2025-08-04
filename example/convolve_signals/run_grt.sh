#!/bin/bash


dist=10
depsrc=2
deprcv=0.5

nt=1024
dt=0.01

modname="milrow"
out="GRN"

grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist}

# convolve different signals
G=$out/${modname}_${depsrc}_${deprcv}_${dist}
P=cout
S=1e24
az=39.2
grt syn -G$G -Osyn_ex -P$P -A$az -S$S -Dt/0.2/0.2/0.4

fn=2
fe=-1
fz=4
grt syn -G$G -Osyn_sf -P$P -A$az -S$S -F$fn/$fe/$fz -Dt/0.1/0.3/0.6

stk=77
dip=88
rak=99
grt syn -G$G -Osyn_dc -P$P -A$az -S$S -M$stk/$dip/$rak -Dp/0.6

M11=1
M12=-2
M13=-5
M22=0.5
M23=3
M33=1.2
grt syn -G$G -Osyn_mt -P$P -A$az -S$S -T$M11/$M12/$M13/$M22/$M23/$M33 -Dr/3

