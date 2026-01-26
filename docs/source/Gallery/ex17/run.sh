#!/bin/bash

set -euo pipefail

rm -rf GRN* *.svg *.tar.gz

cat > halfspace <<EOF
0.0   8.0   4.62    3.3
EOF

cat > halfspace_Q <<EOF
0.0   8.0   4.62    3.3    80    20
EOF

# 强衰减介质需要更高的积分上限
K="-K+k60+s2"
grt greenfn -Mhalfspace   -N2000/0.001+a -D0/0  -R5  -OGRN     $K -Cp
grt greenfn -Mhalfspace_Q -N2000/0.001+a -D0/0  -R5  -OGRN_Q   $K -Cp
python plot.py

cp compare.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf halfspace* GRN* 



