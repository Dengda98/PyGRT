#!/bin/bash

set -euo pipefail

rm -rf *.svg *.tar.gz

python plot.py

cp lamb1_compare.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

