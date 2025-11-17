#!/bin/bash

set -euo pipefail

grt syn -h

grt greenfn -M../milrow -D2/0 -N600/0.02 -R10 -e -OGRN

grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -S1e16 -F-1/2/-4 -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -M33/44/55 -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -Su1e10 -M33/44/55 -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20  -T1/-2/-5/0.5/3/1.2 -Osyn 

grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -Dp/0.6 -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -Dt/0.2/0.4/0.7 -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -Dt/0.4/0.4/0.8 -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -Dr/1.2 -Osyn 
cat > tfile <<EOF
0       0.0
0.02    0.1
0.04    0.2
0.06    0.4
0.08    0.4
0.10    0.4
0.12    0.2
0.14    0.1
0.16    0.0
EOF
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -D0/tfile -Osyn 
rm tfile

grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -M33/44/55 -I1 -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -M33/44/55 -J1 -Osyn 

grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -N -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -e -Osyn 
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -N -e -Osyn 


python -u test_syn.py


rm -rf GRN syn