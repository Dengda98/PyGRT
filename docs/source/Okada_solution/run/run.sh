#!/bin/bash

set -euo pipefail

rm -rf *.nc *.svg

cat > faults.inp <<EOF
  #   X-start    Y-start     X-fin     Y-fin  Kode  shear(m)  reverse(m)  dip  top  bot
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx  xxxxxxxxxx  xxx  xxx  xxx
  1     0.0000     0.0000    20.0000     0.0000 100     0.1000      0.2000  70.0  3.0 10.0
EOF

# BEGIN POINT
grt okada -I6/3.464/2.7 -Su1e12 -Ds5 -Dr0 -M33/44/55 -N -e -X-5/5/0.2 -Y-5/5/0.2 -Ookada_point.nc
# END POINT

# BEGIN FINITE
grt okada -I6/3.464/2.7 -Cfaults.inp -Dr0 -N -e -X-10/10/0.2 -Y-10/10/0.2 -Ookada_finite.nc
# END FINITE

python plot.py

rm -rf *.nc faults.inp
