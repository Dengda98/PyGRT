#!/bin/bash

set -euo pipefail

grt okada -h

grt okada -I6/3.464/2.7 -Su1e12 -Ds50 -Dr0 \
    -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ex.nc
grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Dr0 \
    -M100/20/80 -N -e -X-5/5/0.5 -Y-5/5/0.5 -Ookada_dc.nc
grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Dr0 \
    -M100/20 -N -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ts.nc

python -u test_okada.py

echo "test_okada.sh: all checks passed"

rm -rf okada_ex.nc okada_dc.nc okada_ts.nc
