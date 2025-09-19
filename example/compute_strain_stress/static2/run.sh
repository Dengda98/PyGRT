#!/bin/bash 


mod="halfspace2"
grt static greenfn -M$mod -D20/16 -X3/3/1 -Y-10/10/0.5 -e -Ostgrn1.nc

grt static syn -M0/70/0 -Su1e6 -e -N -Gstgrn1.nc -Osyn_1_stkslip.nc
grt static stress syn_1_stkslip

grt static syn -M0/70/90 -Su1e6 -e -N -Gstgrn1.nc -Osyn_1_dipslip.nc
grt static stress syn_1_dipslip.nc



grt static greenfn -M$mod -D10/14 -X3/3/1 -Y-10/10/0.5 -e -Ostgrn2.nc

grt static syn -M0/70/0 -Su1e6 -e -N -Gstgrn2.nc -Ostsyn_2_stkslip.nc
grt static stress stsyn_2_stkslip.nc

grt static syn -M0/70/90 -Su1e6 -e -N -Gstgrn2.nc -Ostsyn_2_dipslip.nc
grt static stress stsyn_2_dipslip.nc


mod="mod"
grt static greenfn -M$mod -D5/0 -X3/3/1 -Y-10/10/0.5 -e -Ostgrn3.nc

grt static syn -M0/70/0 -Su1e6 -e -N -Gstgrn3.nc -Ostsyn_3_stkslip.nc
grt static stress stsyn_3_stkslip.nc

grt static syn -M0/70/90 -Su1e6 -e -N -Gstgrn3.nc -Ostsyn_3_dipslip.nc
grt static stress stsyn_3_dipslip.nc
