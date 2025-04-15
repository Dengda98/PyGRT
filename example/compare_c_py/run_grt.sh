#!/bin/bash


dist=10
depsrc=2
deprcv=3.3

nt=1024
dt=0.01

modname="milrow"
out="GRN"


#-------------------------- Dynamic -----------------------------------------
grt -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -e

# convolve different signals
G=$out/${modname}_${depsrc}_${deprcv}_${dist}
P=cout
S=1e24
az=39.2
for N in "-N" "" ; do
grt.syn -G$G -Osyn_exp -P$P -A$az -S$S -Dt/0.2/0.2/0.4 -e $N
grt.strain syn_exp/$P
grt.stress syn_exp/$P

fn=2
fe=-1
fz=4
grt.syn -G$G -Osyn_sf -P$P -A$az -S$S -F$fn/$fe/$fz -Dt/0.1/0.3/0.6 -e $N
grt.strain syn_sf/$P
grt.stress syn_sf/$P

stk=77
dip=88
rak=99
grt.syn -G$G -Osyn_dc -P$P -A$az -S$S -M$stk/$dip/$rak -Dp/0.6 -e $N
grt.strain syn_dc/$P
grt.stress syn_dc/$P

M11=1
M12=-2
M13=-5
M22=0.5
M23=3
M33=1.2
grt.syn -G$G -Osyn_mt -P$P -A$az -S$S -T$M11/$M12/$M13/$M22/$M23/$M33 -Dr/3 -e $N
grt.strain syn_mt/$P
grt.stress syn_mt/$P
done


#-------------------------- Static -----------------------------------------
x1=-3
x2=3
nx=11

y1=-4
y2=4
ny=16

mkdir -p static
cd static
stgrt -M../${modname} -D${depsrc}/${deprcv} -X$x1/$x2/$nx -Y$y1/$y2/$ny -e > grn

for N in "-N" "" ; do
stgrt.syn -S$S -e $N < grn > stsyn_exp$N
stgrt.strain < stsyn_exp$N > strain_exp$N
stgrt.stress < stsyn_exp$N > stress_exp$N

stgrt.syn -S$S -F$fn/$fe/$fz -e $N < grn > stsyn_sf$N
stgrt.strain < stsyn_sf$N > strain_sf$N
stgrt.stress < stsyn_sf$N > stress_sf$N

stgrt.syn -S$S -M$stk/$dip/$rak -e $N < grn > stsyn_dc$N
stgrt.strain < stsyn_dc$N > strain_dc$N
stgrt.stress < stsyn_dc$N > stress_dc$N

stgrt.syn -S$S -T$M11/$M12/$M13/$M22/$M23/$M33 -e $N < grn > stsyn_mt$N
stgrt.strain < stsyn_mt$N > strain_mt$N
stgrt.stress < stsyn_mt$N > stress_mt$N
done

cd -

for p in $(ls static/*); do
    echo "----------------------------- $p -------------------------"
    cat $p
done