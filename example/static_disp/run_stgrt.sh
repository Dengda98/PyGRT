#!/bin/bash 

depsrc=2
deprcv=0

x1=-3
x2=3
nx=20

y1=-3
y2=3
ny=20

stgrt -Mhalfspace2 -D${depsrc}/${deprcv} -X$x1/$x2/$nx -Y$y1/$y2/$ny > grn 

# Fault
S="1e23"
stk=60
dip=90
rak=0
stgrt.syn -S$S -M$stk/$dip/$rak -N < grn > syn 

gmt set FONT_TITLE 9p
gmt begin disp_dc png E300
    gmt xyz2grd syn -GsynZ.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,2 -:
    gmt xyz2grd syn -GsynN.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,3 -:
    gmt xyz2grd syn -GsynE.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,4 -:

    gmt basemap -Baf -BWSen+t"Fault, $stk/$dip/$rak, $S" -JX5c/5c -R$y1/$y2/$x1/$x2
    gmt grdimage synZ.nc
    gmt grdvector synE.nc synN.nc -Q0.1c+e+jc+h1+gblack -Si0.3c

    gmt meca -Sa0.5c <<EOF
0 0 $depsrc $stk $dip $rak 5
EOF

    gmt colorbar -Bx+l"Z (cm)"
gmt end


# Single Force
S="1e17"
fn=1
fe=1
fz=0
stgrt.syn -S$S -F$fn/$fe/$fz -N < grn > syn 
gmt begin disp_sf  png E300
    gmt xyz2grd syn -GsynZ.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,2 -:
    gmt xyz2grd syn -GsynN.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,3 -:
    gmt xyz2grd syn -GsynE.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,4 -:

    gmt basemap -Baf -BWSen+t"Force, $fn/$fe/$fz, $S" -JX5c/5c -R$y1/$y2/$x1/$x2
    gmt grdimage synZ.nc
    gmt grdvector synE.nc synN.nc -Q0.1c+e+jc+h1+gblack -Si0.5c

    gmt colorbar -Bx+l"Z (cm)"
gmt end

# Explosion
S="1e23"
stgrt.syn -S$S -N < grn > syn 
gmt begin disp_exp  png E300
    gmt xyz2grd syn -GsynZ.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,2 -:
    gmt xyz2grd syn -GsynN.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,3 -:
    gmt xyz2grd syn -GsynE.nc -R$y1/$y2/$x1/$x2 -I$ny+n/$nx+n -i0,1,4 -:

    gmt basemap -Baf -BWSen+t"Explosion, $S" -JX5c/5c -R$y1/$y2/$x1/$x2
    gmt grdimage synZ.nc
    gmt grdvector synE.nc synN.nc -Q0.1c+e+jc+h1+gblack -Si0.3c

    gmt colorbar -Bx+l"Z (cm)"
gmt end

