#!/bin/bash 

set -euo pipefail

rm -rf *.nc *.svg *.tar.gz

depsrc=2
deprcv=0

x1=-3
x2=3
dx=0.3

y1=-3
y2=3
dy=0.3

grt static greenfn -Mhalfspace -D${depsrc}/${deprcv} -X$x1/$x2/$dx -Y$y1/$y2/$dy -Ostgrn.nc 

# Fault
S="1e23"
stk=60
dip=90
rak=0
grt static syn -S$S -M$stk/$dip/$rak -N -Gstgrn.nc -Ostsyn.nc 

gmt set FONT_TITLE 9p
gmt begin disp_dc pdf 
    gmt basemap -Baf -BWSen+t"Fault, $stk/$dip/$rak, $S" -JX5c/5c -R$y1/$y2/$x1/$x2
    gmt grdimage stsyn.nc?Z
    gmt grdvector stsyn.nc?E stsyn.nc?N -Q0.1c+e+jc+h1+gblack -Si0.3c

    gmt meca -Sa0.5c <<EOF
0 0 $depsrc $stk $dip $rak 5
EOF

    gmt colorbar -Bx+l"Z (cm)"
gmt end
pdf2svg disp_dc.pdf disp_dc.svg


# Single Force
S="1e17"
fn=1
fe=1
fz=0
grt static syn -S$S -F$fn/$fe/$fz -N -Gstgrn.nc -Ostsyn.nc
gmt begin disp_sf pdf 
    gmt basemap -Baf -BWSen+t"Force, $fn/$fe/$fz, $S" -JX5c/5c -R$y1/$y2/$x1/$x2
    gmt grdimage stsyn.nc?Z
    gmt grdvector stsyn.nc?E stsyn.nc?N -Q0.1c+e+jc+h1+gblack -Si0.5c

    gmt colorbar -Bx+l"Z (cm)"
gmt end
pdf2svg disp_sf.pdf disp_sf.svg

# Explosion
S="1e23"
grt static syn -S$S -N -Gstgrn.nc -Ostsyn.nc 
gmt begin disp_ex pdf 
    gmt basemap -Baf -BWSen+t"Explosion, $S" -JX5c/5c -R$y1/$y2/$x1/$x2
    gmt grdimage stsyn.nc?Z
    gmt grdvector stsyn.nc?E stsyn.nc?N -Q0.1c+e+jc+h1+gblack -Si0.3c

    gmt colorbar -Bx+l"Z (cm)"
gmt end
pdf2svg disp_ex.pdf disp_ex.svg


cp disp_dc.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -


rm -rf *.nc *.pdf gmt.conf


