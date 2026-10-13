#!/bin/bash

set -euo pipefail

grt fault2nc -Ffault.inr -Ofault.nc
# 与 USGS 图相同，将零延迟子断层中心设为断层面坐标原点
gmt grdedit "fault.nc?slip"  -R-145/105/-17.5/47.5 -Gslip.nc
gmt grdedit "fault.nc?delay" -R-145/105/-17.5/47.5 -Gdelay.nc
gmt grd2xyz delay.nc | awk '$3 == 0 { print $1, $2 }' > rupture_origin.txt

# 从全局属性读取完整断层面的顶深和底深
read -r depth_top depth_bot < <(ncdump -h fault.nc | awk '/:top =/ { top = $3 } /:bot =/ { print top, $3 }')

gmt begin slip_delay pdf
    gmt set MAP_FRAME_TYPE plain FONT_ANNOT_PRIMARY 10p FONT_LABEL 11p MAP_TICK_LENGTH_PRIMARY 0.10c
    # 自建连续色表，零滑动量为白色，再经浅黄、橙红渐变为深紫
    gmt makecpt -Cwhite,255/255/205,255/216/152,252/143/97,225/74/103,158/47/127,74/14/112 -T0/8/0.2 -Z -D
    # 按实际分块尺寸填色，两个方向保持相同比例，倾向轴向下
    gmt grdimage slip.nc -R-145/105/-17.5/47.5 -JX18c/-4.68c -C -T+o0.10p,gray75 \
        -Bxa50f10+l"Along strike (km)" -Bya10f5+l"Down dip (km)" -BWSn
    gmt grdcontour delay.nc -C5 -A10+f9p,Helvetica,black+gwhite+u" s" -Gn1 -Wc0.5p,gray35 -Wa0.9p,black
    gmt plot rupture_origin.txt -Sa0.32c -Gcyan -W0.8p,black
    gmt text rupture_origin.txt -F+f9p,Helvetica-Bold,black+jBL+t"0 s" -D0.18c/0.18c -Gwhite@10
    gmt colorbar -C -DJBC+w12c/0.35c+h+o0c/1.15c -Bxa1f0.5+l"Slip (m)"
    # 在相同尺寸的图幅上仅叠加右轴，深度向下增加
    gmt basemap -R-145/105/"$depth_top"/"$depth_bot" -JX18c/-4.68c -Bya5f1+l"Depth (km)" -BE
gmt end

pdf2svg slip_delay.pdf slip_delay.svg

rm -f slip_delay.pdf fault.nc slip.nc delay.nc rupture_origin.txt gmt.history
