#!/bin/bash

set -euo pipefail

rm -f lamb_compare.svg

cat > halfspace <<EOF
0.0   8.0   4.62  3.3
EOF

# BEGIN LAMB
grt lamb -H8/4.62/3.3 -N1300/0.005 -R15 -Ds5 -Dr1 -A30 -S1e24 -M33/50/120 -Dt/0.1/0.1/0.2 -n -Olamb_cli
# END LAMB

python plot_compare.py

rm halfspace
rm -rf lamb_cli
