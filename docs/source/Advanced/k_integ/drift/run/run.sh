#!/bin/bash

set -euo pipefail

rm -rf GRN* pygrtstats *.svg

# 没必要执行，故改为如下
cat > /dev/null <<EOF
# BEGIN 1
grt greenfn -Mmilrow -D10/0 -N500/10+a -R5000 -OGRN -S
# END 1
EOF

python run.py

rm -rf pygrtstats