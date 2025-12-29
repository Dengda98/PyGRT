#!/bin/bash

set -euo pipefail

rm -f *.svg

# -----------------------------------------------------------------------------------
# BEGIN LAMB1
grt lamb1 -P0.25 -T0/2/1e-4 -A30 > lamb1.txt
# END LAMB1
# -----------------------------------------------------------------------------------

head -n 10 lamb1.txt > head_lamb1

python lamb1_plot_time.py
python lamb1_plot_freq_time.py
