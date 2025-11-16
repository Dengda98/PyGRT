#!/bin/bash

set -euo pipefail

grt travt -h 

grt travt -M../milrow -D2/0 -R2,3,4,5
cat > dists <<EOF
2
3
4
5
EOF
grt travt -M../milrow -D2/0 -Rdists
rm dists

python -u test_travt.py