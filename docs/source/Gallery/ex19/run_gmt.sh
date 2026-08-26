#!/bin/bash

set -euo pipefail

function plot_disp(){
    local syn=$1
    local figname=$2

region=$(gmt grdinfo -Ir "$syn?Z")
gmt begin $figname png
    gmt set FONT_TITLE 12p
    gmt set MAP_TITLE_OFFSET 2p
    gmt makecpt -CCoulomb_anatolia.cpt -T-2/2/0.01 -D -Z
    gmt basemap "$region" -JM9c -Baf -BWSen
    gmt grdimage "$syn?Z"
    gmt grdvector "$syn?E" "$syn?N" -S50c+s -Q0.2c+e+jb+h0+gblack+n1c -l"50 cm"
    gmt legend -DjBR+w2.5c+o0.2c -F+gwhite+p0.5p -Mh

    gmt colorbar -Bx+l"Vertical displacement Z (cm)"
gmt end

}


function plot_coulomb(){
    local syn=$1
    local figname=$2

region=$(gmt grdinfo -Ir "$syn?Z")
gmt begin $figname png
    gmt set FONT_TITLE 15p
    gmt set MAP_TITLE_OFFSET 2p
    gmt makecpt -CCoulomb_anatolia.cpt -T-0.3/0.3/0.01 -D -Z

    gmt basemap "$region" -JM7c -Baf -BWSen+t"Shear Stress Change"
    gmt grdmath "$syn?tau_s" 1e-7 MUL = tmp.nc
    gmt grdimage tmp.nc

    gmt basemap "$region" -JM7c -Baf -BWSen+t"Normal Stress Change" -X+w+2c
    gmt grdmath "$syn?sigma_n" 1e-7 MUL = tmp.nc
    gmt grdimage tmp.nc
    gmt colorbar -DJBC+w7c -Bx+l"MPa"

    gmt basemap "$region" -JM7c -Baf -BWSen+t"Coulomb Stress Change" -X+w+2c
    gmt grdmath "$syn?coulomb" 1e-7 MUL = tmp.nc
    gmt grdimage tmp.nc

    rm tmp.nc
gmt end
}



plot_disp stsyn_ff.nc.geo  layer_disp
plot_disp okada_ff.nc.geo  halfspace_disp

plot_coulomb stsyn_ff.nc.geo  layer_coulomb
plot_coulomb okada_ff.nc.geo  halfspace_coulomb



