#!/bin/bash

dist=100
depsrc=20
deprcv=8.2

nt=512
dt=0.2

XF=1.0
XL=$(echo $dist*20 | bc)

modname="liquid"
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
rspec96 

# to ascii
hpulse96 -D -i -IMP > hpulse96.out

mkdir -p GRN
f96tosac -B < hpulse96.out
mv *.sac GRN/
