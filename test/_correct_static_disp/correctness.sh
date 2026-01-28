#!/bin/bash
# Check the correctness

set -euo pipefail

rm -rf *.nc

# static greenfn
grt static greenfn -M../milrow -D2/3 -R0/6/0.2 -Cp -e -Ostgrn.nc
grt static greenfn -M../milrow -D2/3 -X-3.1/3.1/0.2 -Y-2.1/2.1/0.2 -Cp -e -Ostgrn.nc

# syn
grt static syn -S1e20 -e -Gstgrn.nc -Ostsyn_ex.nc
grt static syn -S1e20 -e -F2/-1/4    -Gstgrn.nc -Ostsyn_sf.nc
grt static syn -S1e20 -e -M77/88/111 -Gstgrn.nc -Ostsyn_dc.nc
grt static syn -S1e20 -e -M77/88 -Gstgrn.nc -Ostsyn_ts.nc
grt static syn -S1e20 -e -T1/-2/-5/0.5/3/1.2 -X-3.1/3.1/0.3 -Y-2.1/2.1/0.3 -Gstgrn.nc -Ostsyn_mt.nc # new XY grid
grt static syn -S1e20 -e -T1/-2/-5/0.5/3/1.2 -Gstgrn.nc -Ostsyn_mt.nc

# rotate to ZNE
grt static syn -S1e20 -e -T1/-2/-5/0.5/3/1.2 -N -Gstgrn.nc -Ostsyn_mt_ZNE.nc

# strain, stress, rotation
grt static strain stsyn_dc.nc
grt static stress stsyn_dc.nc
grt static rotation stsyn_dc.nc

grt static strain stsyn_mt_ZNE.nc
grt static stress stsyn_mt_ZNE.nc
grt static rotation stsyn_mt_ZNE.nc



python ../compare_nc.py stgrn.nc _Ref/stgrn.nc
python ../compare_nc.py stsyn_ex.nc _Ref/stsyn_ex.nc
python ../compare_nc.py stsyn_sf.nc _Ref/stsyn_sf.nc
python ../compare_nc.py stsyn_dc.nc _Ref/stsyn_dc.nc
python ../compare_nc.py stsyn_ts.nc _Ref/stsyn_ts.nc
python ../compare_nc.py stsyn_mt.nc _Ref/stsyn_mt.nc
python ../compare_nc.py stsyn_mt_ZNE.nc _Ref/stsyn_mt_ZNE.nc

rm -rf *.nc

