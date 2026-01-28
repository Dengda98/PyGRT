#!/bin/bash

# 再测试静态解合成时传入新网格的结果

set -euo pipefail

depsrc=2
deprcv=3.3

modname="milrow"
out="GRN"

rm -rf $out syn_*


#-------------------------- Static -----------------------------------------
x1=-3.1
x2=3.1
dx=0.6

y1=-4.1
y2=4.1
dy=0.8

S=1e24
fn=2
fe=-1
fz=4

stk=77
dip=88
rak=99

M11=1
M12=-2
M13=-5
M22=0.5
M23=3
M33=1.2

rm -rf static

mkdir -p static
cd static
grt static greenfn -M../../${modname} -D${depsrc}/${deprcv} -R0/10/0.1 -e -Ostgrn.nc

for N in "" "-N" ; do
grt static syn -S$S -e $N -Gstgrn.nc -Ostsyn_ex$N.nc -X$x1/$x2/$dx -Y$y1/$y2/$dy
grt static strain stsyn_ex$N.nc
grt static rotation stsyn_ex$N.nc
grt static stress stsyn_ex$N.nc

grt static syn -S$S -F$fn/$fe/$fz -e $N -Gstgrn.nc -Ostsyn_sf$N.nc -X$x1/$x2/$dx -Y$y1/$y2/$dy
grt static strain stsyn_sf$N.nc
grt static rotation stsyn_sf$N.nc
grt static stress stsyn_sf$N.nc

grt static syn -S$S -M$stk/$dip/$rak -e $N -Gstgrn.nc -Ostsyn_dc$N.nc -X$x1/$x2/$dx -Y$y1/$y2/$dy
grt static strain stsyn_dc$N.nc
grt static rotation stsyn_dc$N.nc
grt static stress stsyn_dc$N.nc

grt static syn -S$S -M$stk/$dip -e $N -Gstgrn.nc -Ostsyn_ts$N.nc -X$x1/$x2/$dx -Y$y1/$y2/$dy
grt static strain stsyn_ts$N.nc
grt static rotation stsyn_ts$N.nc
grt static stress stsyn_ts$N.nc

grt static syn -S$S -T$M11/$M12/$M13/$M22/$M23/$M33 -e $N -Gstgrn.nc -Ostsyn_mt$N.nc -X$x1/$x2/$dx -Y$y1/$y2/$dy
grt static strain stsyn_mt$N.nc
grt static rotation stsyn_mt$N.nc
grt static stress stsyn_mt$N.nc
done

cd -


# run pygrt and compare
python -u compare_staticXY.py


rm -rf $out syn_*
rm -rf static