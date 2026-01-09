#!/bin/bash

set -euo pipefail

rm -rf *.svg *.tar.gz

python plot.py halfspace
python plot_all.py halfspace
python plot_all.py milrow

cp halfspace.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

