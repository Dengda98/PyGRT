#!/bin/bash

set -euo pipefail

rm -rf *.svg *.tar.gz

python plot_raypaths.py
python plot_waveforms.py

cp ray_paths.svg cover.svg

ex=$(basename "$(pwd)")
cd .. && tar -czvf "${ex}.tar.gz" "${ex}" && mv "${ex}.tar.gz" "${ex}" && cd -
