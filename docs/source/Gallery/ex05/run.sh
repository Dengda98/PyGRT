#!/bin/bash

set -euo pipefail

rm -rf *.tar.gz *.svg

python run.py

cp trig.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -
