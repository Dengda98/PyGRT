#!/bin/bash

set -euo pipefail

rm -rf *.svg *.tar.gz

python run.py

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -


# ==================== grt 程序 ==========================
# dist=180
# depsrc=5
# deprcv=0

# nt=1400
# dt=0.1

# modname="milrow"
# out="GRN"

# grt greenfn -M${modname} -O${out} -N${nt}/${dt} -D${depsrc}/${deprcv} -R${dist}

# G=$out/${modname}_${depsrc}_${deprcv}_${dist}
# S=1e24
# az=39.2

# stk=77
# dip=88
# rak=99
# grt syn -G$G -Osyn_dc -A$az -S$S -M$stk/$dip/$rak

# rm -rf GRN* syn*


