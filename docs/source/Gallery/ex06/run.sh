#!/bin/bash

set -euo pipefail

rm -rf GRN* *.svg *.tar.gz

### DWM take much time
# grt greenfn -Mmilrow -OGRN_DWM -N1400/1 -D2/0 -R1800 -L20 -E100
grt greenfn -Mmilrow -OGRN_SAFIM -N1400/1 -D2/0 -R1800 -L+a1e-2 -E100

python plot.py

mv safim.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf GRN*
