#!/bin/bash

set -euo pipefail

rm -rf GRN* syn* *.svg *.tar.gz

# No Upsampling
grt greenfn -Mmilrow -N200/0.5 -D2/0 -R100 -OGRN

# Upsampling 5 times
grt greenfn -Mmilrow -N200/0.5+n5 -D2/0 -R100 -OGRN_up

python plot.py

cp upsampling.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf GRN* syn*