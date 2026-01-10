#!/bin/bash

set -euo pipefail

rm -rf *.svg *.tar.gz

python plot.py free halfspace n
python plot.py rigid halfspace n
python plot.py free free n
python plot.py halfspace halfspace y
python infinite_space.py

cp free_halfspace.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

