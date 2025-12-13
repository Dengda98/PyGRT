#!/bin/bash

set -euo pipefail

grt greenfn -h

grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -e -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02+w0.6+n10 -R10 -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02 -H1/10 -R10 -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02 -L20 -R10 -OGRN

grt greenfn -M../milrow -D2/3 -N600/0.02 -K+k4+s1.2+e1e-3+v1.5 -R10 -OGRN

grt greenfn -M../milrow -D2/0 -N1400/1 -L+a1e-3 -E-2/9 -R2000 -OGRN
grt greenfn -M../milrow -D2/0 -N1400/1+a -L+a1e-3 -E-2/9 -R2000 -OGRN

grt greenfn -M../milrow -D2/3 -N600/0.02 -L20 -R10 -S -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02 -L20 -R10 -S1,10,20 -OGRN

# multi distance
grt greenfn -M../milrow -D2/0 -N600/0.02 -R6,8,10 -OGRN
cat > dists <<EOF
6
8
10
EOF
grt greenfn -M../milrow -D2/0 -N600/0.02 -Rdists -OGRN
rm dists

python -u test_greenfn.py

rm -rf GRN
rm -rf GRN_grtstats