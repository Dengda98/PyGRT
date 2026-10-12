#!/bin/bash
set -euo pipefail

function gmtplot_static(){
    local syn=$1
    local S=$2
    local title=${3:-" "}
    local region

    region=$(gmt grdinfo -Ir "$syn")

    gmt basemap -Baf -BWSen+t"$title" -JX5c/6c "$region"
    gmt grdimage ${syn}?Z
    gmt grdvector ${syn}?E ${syn}?N -Q0.08c+e+jb+h1+gblack $S
}

gmt begin syn_grid pdf
    gmt set FONT_TITLE 12p
    gmt set MAP_TITLE_OFFSET 2p
    gmt subplot begin 2x4 -Fs5c/6c -M0.2c/1c
        gmt subplot set 0
        gmtplot_static stsyn_ex.nc -Si0.03c
        gmt colorbar -Bx+l"Z (cm)"

        gmt subplot set 1
        gmtplot_static stsyn_sf.nc -Si6.5c "-F1/-0.5/2"
        gmt colorbar -Bx+l"Z (cm)"

        gmt subplot set 2
        gmtplot_static stsyn_dc.nc -Si0.03c "-M33/50/120"
        gmt meca -Sa0.5c <<EOF
0 0 2 33 50 120 5
EOF
        gmt colorbar -Bx+l"Z (cm)"

        gmt subplot set 3
        gmtplot_static stsyn_dc2.nc -Si0.03c "-M33/90/0"
        gmt meca -Sa0.5c <<EOF
0 0 2 33 90 0 5
EOF
        gmt colorbar -Bx+l"Z (cm)"

        # 这里张裂源的沙滩球仅绘制 DC+CLVD 分量

        gmt subplot set 4
        gmtplot_static stsyn_ts.nc -Si0.03c "-M33/50"
        gmt meca -Sz0.5c <<EOF
0 0 2 $(python tension2mt.py 33 50 stgrn.nc) 24
EOF
        gmt colorbar -Bx+l"Z (cm)"

        gmt subplot set 5
        gmtplot_static stsyn_ts2.nc -Si0.03c "-M33/90"
        gmt meca -Sz0.5c <<EOF
0 0 2 $(python tension2mt.py 33 90 stgrn.nc) 24
EOF
        gmt colorbar -Bx+l"Z (cm)"

        gmt subplot set 6
        gmtplot_static stsyn_mt.nc -Si0.02c "-T0.1/-0.2/1.0/0.3/-0.5/-2.0"
        gmt meca -Sm0.5c <<EOF
0 0 2 -2.0 0.1  0.3  1.0  0.5  0.2 24
EOF

# X Y depth mrr   mtt   mff   mrt   mrf   mtf   exp
#           mzz   mxx   myy   mxz  -myz  -mxy

        gmt colorbar -Bx+l"Z (cm)"

        gmt subplot set 7
        gmtplot_static stsyn_mt2.nc -Si0.13c "-T0/-0.2/0/0/0/0"
        gmt meca -Sm0.5c <<EOF
0 0 2 0.0 0.0  0.0  0.0  0.0  0.2 24
EOF
        gmt colorbar -Bx+l"Z (cm)"
    gmt subplot end
gmt end

pdf2svg syn_grid.pdf syn_grid.svg
