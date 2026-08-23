#!/bin/bash
set -euo pipefail

rm -f stgrn.nc stsyn_ff.nc finite_fault.pdf finite_fault.svg fault.svg

# 绘制有限断层
python plot_fault.py

# ------------------------------------------------
# BEGIN GRN
# 设置 -Ds 考虑多个震源深度
# 设置 -Dr0 仅考虑接收点在地表
# 设置 -R 考虑多个震中距
# 设置 -e 增加位移偏导数的计算
grt static greenfn -Mmilrow -Ds0/8/0.5 -Dr0 -R0/60/0.5 -e -Ostgrn.nc
# END GRN
# ------------------------------------------------

# ------------------------------------------------
# BEGIN SYN
# 使用 -C 读取 Coulomb 格式的有限断层
# 使用 -X/-Y 计算二维接收点网格，并使用 -N 输出 ZNE 分量
grt static syn -Gstgrn.nc -Cfault.inp -N -X-20/20/1 -Y-20/20/1 -Ostsyn_ff.nc -e
# END SYN
# ------------------------------------------------

# ------------------------------------------------
# BEGIN GMT
syn=stsyn_ff.nc
region=$(gmt grdinfo -Ir "$syn?Z")

gmt begin disp pdf
    gmt set FONT_TITLE 12p
    gmt set MAP_TITLE_OFFSET 2p
    gmt makecpt -Cjet -T-4/4/0.1 -D -Z
    gmt basemap "$region" -JX12c/12c -Baf -BWSen+t"Finite-fault static displacement"
    gmt grdimage "$syn?Z" -E
    gmt grdvector "$syn?E" "$syn?N" -S80c -Q0.2c+e+jb+h1+gblack

    gmt plot -W4p <<EOF
-10 -6
10 6
EOF

    gmt colorbar -Bx+l"Vertical displacement Z (cm)"
gmt end

pdf2svg disp.pdf disp.svg
rm -f disp.pdf
# END GMT
# ------------------------------------------------

# 删除中间 NetCDF，仅保留用于文档的 SVG 图
rm -f stgrn.nc stsyn_ff.nc
