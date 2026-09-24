#!/bin/bash

set -euo pipefail

rm -rf *.svg

# -----------------------------------------------------------------------------------
# BEGIN LAMB1
grt lamb1 -P0.25 -T0/2/1e-4 -A30 > lamb1.txt
# END LAMB1
# -----------------------------------------------------------------------------------

# BEGIN LAMB1 MOVING
grt lamb1 -P0.25 -T0/2/1e-3 -A30 -C0.3 > lamb1_moving.txt
# END LAMB1 MOVING

head -n 1 lamb1.txt > head_lamb1
head -n 1 lamb1_moving.txt > head_lamb1_moving

python lamb1_plot_time.py
python lamb1_plot_moving_source_9_4_1.py
python lamb1_plot_moving_source_9_4_4.py
python lamb1_plot_moving_source_9_4_5.py
python lamb1_plot_freq_time.py
