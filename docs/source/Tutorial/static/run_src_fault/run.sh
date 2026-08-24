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


# ------------------------------------------------
# BEGIN COULOMB
# 计算应力张量
grt static stress stsyn_ff.nc
# 将应力张量投影到指定形态的断层面上
grt static sproj -Gstsyn_ff.nc -M59/90/180
# 指定等效摩擦系数，计算库伦应力
grt static coulomb -Gstsyn_ff.nc -F0.75
# END COULOMB
# ------------------------------------------------

# ------------------------------------------------
# BEGIN GMT COULOMB
syn=stsyn_ff.nc
region=$(gmt grdinfo -Ir "$syn?Z")

gmt begin coulomb pdf
    gmt set FONT_TITLE 15p
    gmt set MAP_TITLE_OFFSET 2p
    gmt makecpt -Cjet -T-0.4/0.4/0.01 -D -Z
    cat > faultline <<EOF
-10 -6
10 6
EOF

    gmt basemap "$region" -JX9c/9c -Baf -BWSen+t"Shear Stress Change"
    gmt grdmath "$syn?tau_s" 1e-7 MUL = tmp.nc
    gmt grdimage tmp.nc -E
    gmt plot -W4p faultline

    gmt basemap "$region" -JX9c/9c -Baf -BWSen+t"Normal Stress Change" -X+w+2c
    gmt grdmath "$syn?sigma_n" 1e-7 MUL = tmp.nc
    gmt grdimage tmp.nc -E
    gmt plot -W4p faultline
    gmt colorbar -DJBC+w9c -Bx+l"MPa"

    gmt basemap "$region" -JX9c/9c -Baf -BWSen+t"Coulomb Stress Change" -X+w+2c
    gmt grdmath "$syn?coulomb" 1e-7 MUL = tmp.nc
    gmt grdimage tmp.nc -E
    gmt plot -W4p faultline

    rm faultline tmp.nc
gmt end

pdf2svg coulomb.pdf coulomb.svg
rm -f coulomb.pdf
# END GMT COULOMB
# ------------------------------------------------




# 删除中间 NetCDF，仅保留用于文档的 SVG 图
rm -f stgrn.nc stsyn_ff.nc
