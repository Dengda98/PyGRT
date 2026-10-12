#!/bin/bash

set -euo pipefail

rm -rf CLI_* PY_* *.svg

# -----------------------------------------------------------------------------------
# BEGIN CLI RCVFN
# P 波使用水平射线参数，SV 波使用底部半空间的入射角
grt rcvfn -Mmod1 -P0.03 -TP -N500/0.1 -A5.0 -E-10 -W -OCLI_P
grt rcvfn -Mmod1 -I10/0 -TS -N500/0.1 -A5.0 -E-10 -W -OCLI_S
# END CLI RCVFN
# -----------------------------------------------------------------------------------

python run.py

rm -rf CLI_* PY_*
