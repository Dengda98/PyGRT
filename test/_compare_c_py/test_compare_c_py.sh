#!/bin/bash
# Compare results from C API and Python API

set -euo pipefail

dist=10
depsrc=2
deprcv=3.3

nt=1024
dt=0.01

modname="milrow"
out="GRN"

rm -rf $out syn_*


#-------------------------- Dynamic -----------------------------------------
grt greenfn -M../${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist} -e

# convolve different signals
G=$out/${modname}_${depsrc}_${deprcv}_${dist}
S=1e24
az=39.2
for N in "" "-N" ; do
grt syn -G$G -Osyn_ex$N -A$az -S$S -Dt/0.2/0.2/0.4 -e $N
grt strain syn_ex$N
grt rotation syn_ex$N
grt stress syn_ex$N

fn=2
fe=-1
fz=4
grt syn -G$G -Osyn_sf$N -A$az -S$S -F$fn/$fe/$fz -Dt/0.1/0.3/0.6 -e $N
grt strain syn_sf$N
grt rotation syn_sf$N
grt stress syn_sf$N

stk=77
dip=88
rak=99
grt syn -G$G -Osyn_dc$N -A$az -S$S -M$stk/$dip/$rak -Dp/0.6 -e $N
grt strain syn_dc$N
grt rotation syn_dc$N
grt stress syn_dc$N

grt syn -G$G -Osyn_ts$N -A$az -S$S -M$stk/$dip -Dp/0.6 -e $N
grt strain syn_ts$N
grt rotation syn_ts$N
grt stress syn_ts$N

M11=1
M12=-2
M13=-5
M22=0.5
M23=3
M33=1.2
grt syn -G$G -Osyn_mt$N -A$az -S$S -T$M11/$M12/$M13/$M22/$M23/$M33 -Dr/3 -e $N
grt strain syn_mt$N
grt rotation syn_mt$N
grt stress syn_mt$N
done


#-------------------------- Static -----------------------------------------
x1=-3.1
x2=3.1
dx=0.6

y1=-4.1
y2=4.1
dy=0.8

rm -rf static

mkdir -p static
cd static
grt static greenfn -M../../${modname} -D${depsrc}/${deprcv} -X$x1/$x2/$dx -Y$y1/$y2/$dy -e -Ostgrn.nc

for N in "" "-N" ; do
grt static syn -S$S -e $N -Gstgrn.nc -Ostsyn_ex$N.nc
grt static strain stsyn_ex$N.nc
grt static rotation stsyn_ex$N.nc
grt static stress stsyn_ex$N.nc

grt static syn -S$S -F$fn/$fe/$fz -e $N -Gstgrn.nc -Ostsyn_sf$N.nc
grt static strain stsyn_sf$N.nc
grt static rotation stsyn_sf$N.nc
grt static stress stsyn_sf$N.nc

grt static syn -S$S -M$stk/$dip/$rak -e $N -Gstgrn.nc -Ostsyn_dc$N.nc
grt static strain stsyn_dc$N.nc
grt static rotation stsyn_dc$N.nc
grt static stress stsyn_dc$N.nc

grt static syn -S$S -M$stk/$dip -e $N -Gstgrn.nc -Ostsyn_ts$N.nc
grt static strain stsyn_ts$N.nc
grt static rotation stsyn_ts$N.nc
grt static stress stsyn_ts$N.nc

grt static syn -S$S -T$M11/$M12/$M13/$M22/$M23/$M33 -e $N -Gstgrn.nc -Ostsyn_mt$N.nc
grt static strain stsyn_mt$N.nc
grt static rotation stsyn_mt$N.nc
grt static stress stsyn_mt$N.nc
done

cd -


# run pygrt and compare
python -u compare.py


rm -rf $out syn_*
rm -rf static