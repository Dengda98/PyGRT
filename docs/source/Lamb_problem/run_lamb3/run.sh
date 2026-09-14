#!/bin/bash

set -euo pipefail

rm -rf *.svg

# BEGIN LAMB3
# grt lamb3 -P0.25 -T0/2/5e-3 -R10 -D2/1 -A30 > lamb3.txt
# 还可以使用 -S 来输出源点坐标偏导、接收点坐标偏导和二阶混合偏导
grt lamb3 -P0.25 -T0/2/5e-3 -R10 -D2/1 -A30 -S+slamb3_source.txt+rlamb3_receiver.txt+mlamb3_mixed.txt > lamb3.txt
# END LAMB3

head -n 1 lamb3.txt > head_lamb3
head -n 1 lamb3_source.txt > head_lamb3_source
head -n 1 lamb3_receiver.txt > head_lamb3_receiver
head -n 1 lamb3_mixed.txt > head_lamb3_mixed

python lamb3_plot_time.py
python lamb3_plot_freq_time.py
python lamb3_plot_freq_time_source.py
python lamb3_plot_freq_time_mixed.py
