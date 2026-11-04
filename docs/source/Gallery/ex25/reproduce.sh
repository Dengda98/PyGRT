#!/bin/bash

set -euo pipefail

python 01_create_input_files.py
python 02_simulate_ambient_noise.py
python 03_segment_ambient_noise.py
python 04_compute_ccf.py
python 05_plot_geom_disp.py
python 06_plot_ccf.py

cp ambient_figures/profiles/STA001.png cover.png
