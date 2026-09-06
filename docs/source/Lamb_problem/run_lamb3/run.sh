#!/bin/bash

set -euo pipefail

rm -rf *.svg

# BEGIN LAMB3
grt lamb3 -P0.25 -T0/2/5e-3 -R10 -D2/1 -A30 -S+slamb3_source.txt+rlamb3_receiver.txt > lamb3.txt
# END LAMB3

head -n 10 lamb3.txt > head_lamb3
head -n 10 lamb3_source.txt > head_lamb3_source
head -n 10 lamb3_receiver.txt > head_lamb3_receiver
echo "..." >> head_lamb3
echo "..." >> head_lamb3_source
echo "..." >> head_lamb3_receiver

python lamb3_plot_time.py
python lamb3_plot_freq_time.py
