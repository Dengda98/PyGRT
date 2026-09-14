#!/bin/bash

set -euo pipefail

rm -rf *.svg
rm -f lamb2_source.txt lamb2_receiver.txt lamb2_mixed.txt

# BEGIN LAMB2
# grt lamb2 -P0.25 -T0/2/1e-3 -R10 -Ds5 -A30 > lamb2.txt
# 还可以使用 -S 来输出源点坐标偏导、接收点坐标偏导和二阶混合偏导
grt lamb2 -P0.25 -T0/2/1e-3 -R10 -Ds5 -A30 -S+slamb2_source.txt+rlamb2_receiver.txt+mlamb2_mixed.txt > lamb2.txt
# END LAMB2

head -n 1 lamb2.txt > head_lamb2
head -n 1 lamb2_source.txt > head_lamb2_source
head -n 1 lamb2_receiver.txt > head_lamb2_receiver
head -n 1 lamb2_mixed.txt > head_lamb2_mixed

python lamb2_plot_time.py
python lamb2_plot_freq_time.py
python lamb2_plot_freq_time_source.py
python lamb2_plot_freq_time_mixed.py
