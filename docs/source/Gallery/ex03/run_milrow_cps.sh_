#!/bin/bash


dist=10
depsrc=2
deprcv=0

nt=1024
dt=0.01

XF=1.0
XL=$(echo $dist*20 | bc)

modname="milrow"
mod="${modname}.mod"

maindir=${modname}_sdep${depsrc}_rdep${deprcv}
mkdir -p ${maindir}
cd ${maindir}
# mkdir -p ${modname}
# cd ${modname}

cat > dfile << EOF
$dist $dt $nt 0 0
EOF


# compute Green's Functions
hprep96 -M ../$mod -d dfile -HS $depsrc -HR $deprcv -ALL -XL $XL -XF $XF
hspec96 

# to ascii
hpulse96 -D -i -IMP > hpulse96.out

mkdir -p GRN
f96tosac -B < hpulse96.out
mv *.sac GRN/
