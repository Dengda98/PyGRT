#!/bin/bash
set -euo pipefail

rm -rf stgrn* stsyn* *.svg

function _example_gmtplot(){
# ---------------------------------------------------------------------------------
# BEGIN gmt
syn="stsyn_ex.nc"
gmt begin syn_ex pdf
    gmt basemap -Baf -BWSen -JX5c/6c $(gmt grdinfo -Ir $syn)
    gmt grdimage ${syn}?Z
    gmt grdvector ${syn}?E ${syn}?N -Q0.08c+e+jc+h1+gblack -Si0.03c
    gmt colorbar -Bx+l"Z (cm)"
gmt end
# END gmt
# ---------------------------------------------------------------------------------
}


function gmtplot_static(){
    local syn=$1
    local S=$2

    gmt basemap -Baf -BWSen -JX5c/6c $(gmt grdinfo -Ir $syn)
    gmt grdimage ${syn}?Z
    gmt grdvector ${syn}?E ${syn}?N -Q0.08c+e+jc+h1+gblack $S
}


# ---------------------------------------------------------------------------------
# BEGIN GRN
depsrc=2
deprcv=0

x1=-3
x2=3
dx=0.15

y1=-2.5
y2=2.5
dy=0.15
# 输出到NetCDF网格文件
grt static greenfn -Mmilrow -D${depsrc}/${deprcv} -X$x1/$x2/$dx -Y$y1/$y2/$dy -Ostgrn.nc
# END GRN
# ---------------------------------------------------------------------------------

ncdump -h stgrn.nc > grn_head


# ---------------------------------------------------------------------------------
# BEGIN SYN EXP
# 从网格文件中读取格林函数，再将合成结果写入新网格
grt static syn -S1e24 -N -Gstgrn.nc -Ostsyn_ex.nc
# END SYN EXP
# ---------------------------------------------------------------------------------

gmt begin syn_ex pdf
    gmtplot_static stsyn_ex.nc -Si0.03c
    gmt colorbar -Bx+l"Z (cm)"
gmt end


# ---------------------------------------------------------------------------------
# BEGIN SYN SF
# 从网格文件中读取格林函数，再将合成结果写入新网格
grt static syn -S1e16 -F1/-0.5/2 -N -Gstgrn.nc -Ostsyn_sf.nc
# END SYN SF
# ---------------------------------------------------------------------------------

gmt begin syn_sf pdf
    gmtplot_static stsyn_sf.nc -Si6.5c
    gmt colorbar -Bx+l"Z (cm)"
gmt end

# ---------------------------------------------------------------------------------
# BEGIN SYN DC
# 从网格文件中读取格林函数，再将合成结果写入新网格
grt static syn -S1e24 -M33/50/120 -N -Gstgrn.nc -Ostsyn_dc.nc
# END SYN DC
# ---------------------------------------------------------------------------------

gmt begin syn_dc pdf
    gmtplot_static stsyn_dc.nc -Si0.03c
    gmt meca -Sa0.5c <<EOF
0 0 2 33 50 120 5
EOF
    gmt colorbar -Bx+l"Z (cm)"
gmt end

# ---------------------------------------------------------------------------------
# BEGIN SYN DC2
# 从网格文件中读取格林函数，再将合成结果写入新网格
grt static syn -S1e24 -M33/90/0 -N -Gstgrn.nc -Ostsyn_dc2.nc
# END SYN DC2
# ---------------------------------------------------------------------------------


gmt begin syn_dc2 pdf
    gmtplot_static stsyn_dc2.nc -Si0.03c
    gmt meca -Sa0.5c <<EOF
0 0 2 33 90 0 5
EOF
    gmt colorbar -Bx+l"Z (cm)"
gmt end

# ---------------------------------------------------------------------------------
# BEGIN SYN MT
# 从网格文件中读取格林函数，再将合成结果写入新网格
grt static syn -S1e24 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -N -Gstgrn.nc -Ostsyn_mt.nc
# END SYN MT
# ---------------------------------------------------------------------------------

gmt begin syn_mt pdf
    gmtplot_static stsyn_mt.nc -Si0.02c
    gmt meca -Sm0.5c <<EOF
0 0 2 -2.0 0.1  0.3  1.0  0.5  0.2 24
EOF
    gmt colorbar -Bx+l"Z (cm)"
gmt end


# ---------------------------------------------------------------------------------
# BEGIN SYN MT2
# 从网格文件中读取格林函数，再将合成结果写入新网格
grt static syn -S1e24 -T0/-0.2/0/0/0/0 -N -Gstgrn.nc -Ostsyn_mt2.nc
# END SYN MT2
# ---------------------------------------------------------------------------------


# X Y depth mrr   mtt   mff   mrt   mrf   mtf   exp
#           mzz   mxx   myy   mxz  -myz  -mxy
gmt begin syn_mt2 pdf
    gmtplot_static stsyn_mt2.nc -Si0.13c
    gmt meca -Sm0.5c <<EOF
0 0 2 0.0 0.0  0.0  0.0  0.0  0.2 24
EOF
    gmt colorbar -Bx+l"Z (cm)"
gmt end

# 统一转为 svg
for pdfname in $(ls *.pdf); do
    name=$(basename $pdfname .pdf)
    pdf2svg $pdfname ${name}.svg
    rm -r $pdfname
done

rm -rf stgrn* stsyn*