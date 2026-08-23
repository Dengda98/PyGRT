#!/bin/bash

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

cat > cfaults.inp <<'EOF'
  #   X-start    Y-start     X-fin     Y-fin    Kode  shear(m)  reverse(m)  dip angle   top(km)   bot(km)
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000     0.0000    30.0000     0.0000 100     0.1000      0.2000  70.0 20.0 30.0
EOF

cat > rcv_pts_6.txt <<'EOF'
# north east depth strike dip rake
0 0 0 10 20 30
1 2 0 40 50 60
-1 1 0 70 80 90
EOF

cat > rcv_faults.inp <<'EOF'
  #   X-start    Y-start      X-fin      Y-fin   Kode  value1      value2       dip       top       bot
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000     0.0000     2.0000     0.0000   100     0.1000      0.0000     90.00       1.0000     3.0000
  2     0.0000     4.0000     2.0000     4.0000   200     0.0000      0.1000     60.00       1.0000     3.0000
EOF

grt okada -h

grt okada -I6/3.464/2.7 -Su1e12 -Ds50 -Dr0 -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ex.nc
grt okada -I6/3.464/2.7 -Su1e12 -Ds10 -Dr10 -X0/0/1 -Y0/0/1 -Ookada_singular.nc
grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Dr0 -M100/20/80 -N -e -X-5/5/0.5 -Y-5/5/0.5 -Ookada_dc.nc
grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Dr0 -M100/20 -N -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ts.nc
grt okada -I6/3.464/2.7 -Su1e16 -Ds10 -Qrcv_pts_6.txt -N -Ookada_q6.nc
grt okada -I6/3.464/2.7 -Ccfaults.inp -Dr0 -e -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ff_zrt_cli.nc
grt okada -I6/3.464/2.7 -Ccfaults.inp -Dr0 -N -e -X-5/5/0.5 -Y-5/5/0.5 -Ookada_ff_zne_cli.nc
grt okada -I6/3.464/2.7 -Su1e12 -Ds10 -Rrcv_faults.inp+i0.75/0.75 -N -e -Ookada_rf.nc
grt okada -I6/3.464/2.7 -Su1e12 -Ds10 -Rrcv_faults.inp -N -e -Ookada_rf_default.nc

python -u test_okada.py

echo "test_okada.sh: all checks passed"

rm -rf okada_ex.nc okada_singular.nc okada_dc.nc okada_ts.nc okada_ff_zrt_cli.nc okada_ff_zne_cli.nc \
    okada_q6.nc okada_rf.nc okada_rf_default.nc okada_python.nc okada_python_ff.nc okada_python_ff_zrt.nc okada_python_ff_zne.nc \
    cfaults.inp rcv_pts_6.txt rcv_faults.inp
