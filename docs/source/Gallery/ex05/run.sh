#!/bin/bash

set -euo pipefail

rm -rf *.tar.gz *.svg

python run.py

cp trig.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -


# dist=10
# depsrc=2
# deprcv=0.5

# nt=1024
# dt=0.01

# modname="milrow"
# out="GRN"

# grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist}

# # convolve different signals
# G=$out/${modname}_${depsrc}_${deprcv}_${dist}
# S=1e24
# az=30

# fn=2
# fe=-1
# fz=4
# grt syn -G$G -Osyn_sf1 -A$az -S$S -F$fn/$fe/$fz -Dt/0.1/0.3/0.6
# grt syn -G$G -Osyn_sf2 -A$az -S$S -F$fn/$fe/$fz -Dp/0.6
# grt syn -G$G -Osyn_sf3 -A$az -S$S -F$fn/$fe/$fz -Dr/3


