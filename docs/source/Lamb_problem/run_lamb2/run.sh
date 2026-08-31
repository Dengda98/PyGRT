#!/bin/bash

set -euo pipefail

rm -rf *.svg

# BEGIN LAMB2
grt lamb2 -P0.25 -T0/2/1e-3 -D60 -A30 > lamb2.txt
# END LAMB2

head -n 10 lamb2.txt > head_lamb2


python lamb2_plot_time.py
python lamb2_plot_freq_time.py
