#!/bin/bash

set -euo pipefail

rm -rf *.svg
rm -f lamb2_source.txt lamb2_receiver.txt

# BEGIN LAMB2
grt lamb2 -P0.25 -T0/2/1e-3 -R10 -Ds5 -A30 -S+slamb2_source.txt+rlamb2_receiver.txt > lamb2.txt
# END LAMB2

head -n 10 lamb2.txt > head_lamb2
head -n 10 lamb2_source.txt > head_lamb2_source
head -n 10 lamb2_receiver.txt > head_lamb2_receiver
echo "..." >> head_lamb2
echo "..." >> head_lamb2_source
echo "..." >> head_lamb2_receiver


python lamb2_plot_time.py
python lamb2_plot_freq_time.py
