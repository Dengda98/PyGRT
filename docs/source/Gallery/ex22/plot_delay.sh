#!/bin/bash

set -euo pipefail

# 将最细剖分的破裂延迟转换为断层面网格
grt fault2nc -Finputs/strike_slip_dx0.5.inp -Oresults/delay_strike_slip.nc
grt fault2nc -Finputs/thrust_dx0.5.inp      -Oresults/delay_thrust.nc

gmt begin delay png
    gmt set FONT_ANNOT_PRIMARY 11p FONT_LABEL 12p FONT_TITLE 14p MAP_FRAME_TYPE plain
    gmt subplot begin 1x2 -Fs10c/5c -M0.8c -R0/20/0/10 -JX10c/-5c
        # 横轴沿走向，纵轴向下倾，两幅图保持相同的距离比例
        gmt subplot set 0
        gmt basemap -Bxa5f1+l"Along strike (km)" -Bya2f1+l"Down dip (km)" -BWSen+t"Strike-slip fault: delay (s)"
        gmt grdcontour "results/delay_strike_slip.nc?delay" -C0.25 -A0.5+f10p -Gn1 -Wc0.6p,gray40 -Wa1p,black

        # 每 0.25 s 绘制等值线，每 0.5 s 标注延迟时间
        gmt subplot set 1
        gmt basemap -Bxa5f1+l"Along strike (km)" -Bya2f1+l"Down dip (km)" -BWSen+t"Reverse fault: delay (s)"
        gmt grdcontour "results/delay_thrust.nc?delay" -C0.25 -A0.5+f10p -Gn1 -Wc0.6p,gray40 -Wa1p,black
    gmt subplot end
gmt end
