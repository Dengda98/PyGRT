#!/bin/bash
set -eu

rm -rf grn* syn* *.png

# ---------------------------------------------------------------------------------
# BEGIN gmt
function gmtplot_static(){
    local syn=$1
    local S=$2

    gmt xyz2grd $syn -GsynZ.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,2 -:
    gmt xyz2grd $syn -GsynN.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,3 -:
    gmt xyz2grd $syn -GsynE.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,4 -:

    gmt basemap -Baf -BWSen -JX5c/5c -R$y1/$y2/$x1/$x2
    gmt grdimage synZ.nc
    gmt grdvector synE.nc synN.nc -Q0.08c+e+jc+h1+gblack $S

    rm *.nc
}
# END gmt
# ---------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
# BEGIN GRN
depsrc=2
deprcv=0

x1=-3
x2=3
nx=41

y1=-2.5
y2=2.5
ny=33
# 输出到标准输出，重定向到grn文件中
stgrt -Mmilrow -D${depsrc}/${deprcv} -X$x1/$x2/$nx -Y$y1/$y2/$ny > grn
# END GRN
# ---------------------------------------------------------------------------------

head -n6 grn > grn_head
echo "..." >> grn_head


# ---------------------------------------------------------------------------------
# BEGIN SYN EXP
# 从标准输入中读取格林函数
stgrt.syn -S1e24 -N < grn > syn_ex
# END SYN EXP
# ---------------------------------------------------------------------------------

# gmt begin syn_ex png E300
#     gmtplot_static syn_ex -Si0.03c
#     gmt colorbar -Bx+l"Z (cm)"
# gmt end


# ---------------------------------------------------------------------------------
# BEGIN SYN SF
# 从标准输入中读取格林函数
stgrt.syn -S1e16 -F1/-0.5/2 -N < grn > syn_sf
# END SYN SF
# ---------------------------------------------------------------------------------

# gmt begin syn_sf png E300
#     gmtplot_static syn_sf -Si6.5c
#     gmt colorbar -Bx+l"Z (cm)"
# gmt end

# ---------------------------------------------------------------------------------
# BEGIN SYN DC
# 从标准输入中读取格林函数
stgrt.syn -S1e24 -M33/50/120 -N < grn > syn_dc
# END SYN DC
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN SYN DC2
# 从标准输入中读取格林函数
stgrt.syn -S1e24 -M33/90/0 -N < grn > syn_dc2
# END SYN DC2
# ---------------------------------------------------------------------------------


# gmt begin syn_dc png E300
#     gmtplot_static syn_dc -Si0.03c
#     gmt meca -Sa0.5c <<EOF
# 0 0 2 33 50 120 5
# EOF
#     gmt colorbar -Bx+l"Z (cm)"
# gmt end

# ---------------------------------------------------------------------------------
# BEGIN SYN MT
# 从标准输入中读取格林函数
stgrt.syn -S1e24 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -N < grn > syn_mt
# END SYN MT
# ---------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------
# BEGIN SYN MT2
# 从标准输入中读取格林函数
stgrt.syn -S1e24 -T0/-0.2/0/0/0/0 -N < grn > syn_mt2
# END SYN MT2
# ---------------------------------------------------------------------------------


# X Y depth mrr   mtt   mff   mrt   mrf   mtf   exp
#           mzz   mxx   myy   mxz  -myz  -mxy
# gmt begin syn_mt png E300
#     gmtplot_static syn_mt -Si0.02c
#     gmt meca -Sm0.5c <<EOF
# 0 0 2 -2.0 0.1  0.3  1.0  0.5  0.2 24
# EOF
#     gmt colorbar -Bx+l"Z (cm)"
# gmt end

