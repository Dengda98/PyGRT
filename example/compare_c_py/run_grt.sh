#!/bin/bash
set -eu

dist=10
depsrc=2
deprcv=3.3

nt=1024
dt=0.01

modname="milrow"
out="GRN"


#-------------------------- Dynamic -----------------------------------------
grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -e

# convolve different signals
G=$out/${modname}_${depsrc}_${deprcv}_${dist}
P=cout
S=1e24
az=39.2
for N in "" "-N" ; do
grt syn -G$G -Osyn_ex -P$P -A$az -S$S -Dt/0.2/0.2/0.4 -e $N
grt strain syn_ex/$P
grt rotation syn_ex/$P
grt stress syn_ex/$P

fn=2
fe=-1
fz=4
grt syn -G$G -Osyn_sf -P$P -A$az -S$S -F$fn/$fe/$fz -Dt/0.1/0.3/0.6 -e $N
grt strain syn_sf/$P
grt rotation syn_sf/$P
grt stress syn_sf/$P

stk=77
dip=88
rak=99
grt syn -G$G -Osyn_dc -P$P -A$az -S$S -M$stk/$dip/$rak -Dp/0.6 -e $N
grt strain syn_dc/$P
grt rotation syn_dc/$P
grt stress syn_dc/$P

M11=1
M12=-2
M13=-5
M22=0.5
M23=3
M33=1.2
grt syn -G$G -Osyn_mt -P$P -A$az -S$S -T$M11/$M12/$M13/$M22/$M23/$M33 -Dr/3 -e $N
grt strain syn_mt/$P
grt rotation syn_mt/$P
grt stress syn_mt/$P
done


#-------------------------- Static -----------------------------------------
x1=-3
x2=3
dx=0.6

y1=-4
y2=4
dy=0.8

mkdir -p static
cd static
grt static greenfn -M../${modname} -D${depsrc}/${deprcv} -X$x1/$x2/$dx -Y$y1/$y2/$dy -e > grn

for N in "" "-N" ; do
grt static syn -S$S -e $N < grn > stsyn_ex$N
grt static strain < stsyn_ex$N > strain_ex$N
grt static rotation < stsyn_ex$N > rotation_ex$N
grt static stress < stsyn_ex$N > stress_ex$N

grt static syn -S$S -F$fn/$fe/$fz -e $N < grn > stsyn_sf$N
grt static strain < stsyn_sf$N > strain_sf$N
grt static rotation < stsyn_sf$N > rotation_sf$N
grt static stress < stsyn_sf$N > stress_sf$N

grt static syn -S$S -M$stk/$dip/$rak -e $N < grn > stsyn_dc$N
grt static strain < stsyn_dc$N > strain_dc$N
grt static rotation < stsyn_dc$N > rotation_dc$N
grt static stress < stsyn_dc$N > stress_dc$N

grt static syn -S$S -T$M11/$M12/$M13/$M22/$M23/$M33 -e $N < grn > stsyn_mt$N
grt static strain < stsyn_mt$N > strain_mt$N
grt static rotation < stsyn_mt$N > rotation_mt$N
grt static stress < stsyn_mt$N > stress_mt$N
done

cd -

for p in $(ls static/*); do
    echo "----------------------------- $p -------------------------"
    cat $p
done