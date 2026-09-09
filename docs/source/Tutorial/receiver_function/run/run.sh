#!/bin/bash

set -euo pipefail

rm -rf CLI_P CLI_S PY_P PY_S rftn.svg rftn_P.svg rftn_S.svg

# -----------------------------------------------------------------------------------
# BEGIN CLI RFTN
# P 波使用水平射线参数，SV 波使用底部半空间的入射角
grt rftn -Mmod1 -P0.03 -TP -N500/0.1 -A5.0 -E-10 -W -OCLI_P
grt rftn -Mmod1 -I10/0 -TS -N500/0.1 -A5.0 -E-10 -W -OCLI_S
# END CLI RFTN
# -----------------------------------------------------------------------------------

python run.py

rm -rf CLI_* PY_*