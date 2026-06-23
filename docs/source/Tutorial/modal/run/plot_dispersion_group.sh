#!/bin/bash

set -euo pipefail

function _plot(){
    filepath=$1
    title=$2
    args=$3

    # 在 mac 上运行 ncap2 总是需要指定“输出文件”， weird
    JUNK=".JUNK"

    # 最高阶
    maxmode=$(ncap2 -O -s 'print(max(mode))' $filepath $JUNK | awk '{print $NF}')
    # 最高频率
    maxfreq=$(ncap2 -O -s 'print(max(freq))' $filepath $JUNK | awk '{print $NF}')
    # 取合适的速度范围
    minc=$(ncap2 -O -s 'print(min(u))' $filepath $JUNK | awk '{print $NF}')
    maxc=$(ncap2 -O -s 'print(max(u))' $filepath $JUNK | awk '{print $NF}')
    c1=$(echo "scale=1; $minc - ($maxc - $minc) * 0.05" | bc)
    c2=$(echo "scale=1; $maxc + ($maxc - $minc) * 0.05" | bc)

    rm $JUNK

    gmt basemap -R0/$maxfreq/$c1/$c2 -JX10c/7c -Bxa1f+l"Frequency (Hz)" -Bya0.2f+l"Group velocity (km/s)" -BWSen+t"$title" $args
    for i in $(seq 0 $maxmode); do
        grt disp2asc -U$filepath -N$i | gmt plot -i0,3 -W0.5p,blue
    done
}


filename="dispersion_sh_group"
gmt begin $filename pdf
    gmt set FONT_TITLE 12p,0

    _plot group_R.nc Rayleigh ""
    _plot group_L.nc Love -X+w+2c
gmt end

pdf2svg $filename.pdf $filename.svg
rm $filename.pdf