#!/bin/bash 

set -euo pipefail

rm -rf *.nc *.svg *.tar.gz

mod="halfspace"
grt static greenfn -M$mod -D20/16 -Y3/3/1 -X-10/10/0.5 -e -Ostgrn1.nc

grt static syn -M90/70/0 -Su1e6 -e -N -Gstgrn1.nc -Ostsyn_1_stkslip.nc
grt static stress stsyn_1_stkslip.nc

grt static syn -M90/70/90 -Su1e6 -e -N -Gstgrn1.nc -Ostsyn_1_dipslip.nc
grt static stress stsyn_1_dipslip.nc


grt static greenfn -M$mod -D10/14 -Y3/3/1 -X-10/10/0.5 -e -Ostgrn2.nc

grt static syn -M90/70/0 -Su1e6 -e -N -Gstgrn2.nc -Ostsyn_2_stkslip.nc
grt static stress stsyn_2_stkslip.nc

grt static syn -M90/70/90 -Su1e6 -e -N -Gstgrn2.nc -Ostsyn_2_dipslip.nc
grt static stress stsyn_2_dipslip.nc


mod="mod"
grt static greenfn -M$mod -D5/0 -Y3/3/1 -X-10/10/0.5 -e -Ostgrn3.nc

grt static syn -M90/70/0 -Su1e6 -e -N -Gstgrn3.nc -Ostsyn_3_stkslip.nc
grt static stress stsyn_3_stkslip.nc

grt static syn -M90/70/90 -Su1e6 -e -N -Gstgrn3.nc -Ostsyn_3_dipslip.nc
grt static stress stsyn_3_dipslip.nc



python plot.py

cp stress1.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf *.nc