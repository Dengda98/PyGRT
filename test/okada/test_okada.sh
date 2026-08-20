#!/bin/bash

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

cat > cfaults.inp <<'EOF'
  #   X-start    Y-start     X-fin     Y-fin    Kode  shear(m)  reverse(m)  dip angle   top(km)   bot(km)
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000     0.0000    30.0000     0.0000 100     0.1000      0.2000  70.0 20.0 30.0
EOF

grt okada -h

grt okada -I6/3.464/2.7 -Su1e12 -Ds50 -Dr0 -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ex.nc
grt okada -I6/3.464/2.7 -Su1e12 -Ds10 -Dr10 -X0/0/1 -Y0/0/1 -Ookada_singular.nc
grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Dr0 -M100/20/80 -N -e -X-5/5/0.5 -Y-5/5/0.5 -Ookada_dc.nc
grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Dr0 -M100/20 -N -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ts.nc
grt okada -I6/3.464/2.7 -Ccfaults.inp -Dr0 -N -e -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ff.nc

python -u test_okada.py

echo "test_okada.sh: all checks passed"

rm -rf okada_ex.nc okada_singular.nc okada_dc.nc okada_ts.nc okada_ff.nc okada_python.nc okada_python_ff.nc cfaults.inp
