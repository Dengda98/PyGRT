#!/bin/bash

set -euo pipefail

grt strain -h
grt stress -h
grt rotation -h

grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -e -OGRN
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -e -Osyn 
grt strain syn
grt stress syn
grt rotation syn

grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -e -N -Osyn_ZNE 
grt strain syn_ZNE 
grt stress syn_ZNE 
grt rotation syn_ZNE 

# -------------------- 静态应变 / 应力 / 旋转 --------------------
grt static strain -h
grt static stress -h
grt static rotation -h

grt static greenfn -M../milrow -D2/0 -X-3/3/1 -Y-2/2/1 -e -Ostgrn.nc
grt static syn -Gstgrn.nc -S1e20 -e -Ostsyn.nc
grt static strain stsyn.nc
grt static stress stsyn.nc
grt static rotation stsyn.nc

grt static syn -Gstgrn.nc -S1e20 -e -N -Ostsyn_ZNE.nc
grt static strain stsyn_ZNE.nc
grt static stress stsyn_ZNE.nc
grt static rotation stsyn_ZNE.nc

# 任意接收点布局
cat > rcv_pts.txt <<'EOF'
# north east depth (km)
0 0 0
1 2 0
-1 1 0
EOF
grt static syn -Gstgrn.nc -S1e20 -e -Qrcv_pts.txt -Ostsyn_q.nc
grt static strain stsyn_q.nc
grt static stress stsyn_q.nc
grt static rotation stsyn_q.nc

python -u test_tensors.py

rm -rf GRN syn* stgrn.nc stsyn*.nc rcv_pts.txt
