#!/bin/bash
set -euo pipefail

rm -rf stgrn* stsyn*
rm -f syn.svg syn.pdf rcv_pts.txt

# ---------------------------------------------------------------------------------
# BEGIN GRN
grt static greenfn -Mmilrow -D2/0 -R0/10/0.1 -Ostgrn.nc
# END GRN
# ---------------------------------------------------------------------------------

ncdump -h stgrn.nc > grn_head


# ---------------------------------------------------------------------------------
# BEGIN SYN GRID
# 爆炸源
grt static syn -Gstgrn.nc -S1e24 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_ex.nc -N

# 单力源
grt static syn -Gstgrn.nc -S1e16 -F1/-0.5/2 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_sf.nc -N

# 剪切源
grt static syn -Gstgrn.nc -S1e24 -M33/50/120 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_dc.nc -N

# 剪切源
grt static syn -Gstgrn.nc -S1e24 -M33/90/0 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_dc2.nc -N

# 张裂源
grt static syn -Gstgrn.nc -S1e24 -M33/50 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_ts.nc -N

# 张裂源
grt static syn -Gstgrn.nc -S1e24 -M33/90 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_ts2.nc -N

# 矩张量源
grt static syn -Gstgrn.nc -S1e24 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_mt.nc -N

# 矩张量源
grt static syn -Gstgrn.nc -S1e24 -T0/-0.2/0/0/0/0 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_mt2.nc -N
# END SYN GRID
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN GMT
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
# END GMT
# ---------------------------------------------------------------------------------

# 统一转为 svg
for pdfname in $(ls *.pdf); do
    name=$(basename $pdfname .pdf)
    pdf2svg $pdfname ${name}.svg
    rm -rf $pdfname
done


python - <<'PY'
import numpy as np

theta = np.linspace(0.0, 2.0 * np.pi, 72, endpoint=False)
radius = 5.0
points = np.column_stack((
    radius * np.cos(theta),
    radius * np.sin(theta),
    np.zeros_like(theta),
))
np.savetxt("rcv_pts.txt", points, fmt="%.8f", header="north east depth (km)")
PY

head -n 10 rcv_pts.txt > rcv_head
echo "..." >> rcv_head


# ---------------------------------------------------------------------------------
# BEGIN SYN POINTS
# 使用 -Q 来传入任意点坐标文件
grt static syn -Gstgrn.nc -S1e24 -M33/90/0 -Qrcv_pts.txt -Ostsyn_points.nc -N
# END SYN POINTS
# ---------------------------------------------------------------------------------

python plot_points.py

rm -rf stgrn* stsyn* rcv_pts.txt
