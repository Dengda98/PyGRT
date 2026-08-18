#!/bin/bash

set -euo pipefail

rm -rf *.nc *.svg *.tar.gz

cat > halfspace <<EOF
0.0   6.0   3.464   2.7
EOF

depsrc=5
deprcv=4
x1=-5
x2=5
dx=0.2
y1=-5
y2=5
dy=0.2
strike=33
dip=44
rake=55

grt static greenfn -Mhalfspace -D$depsrc/$deprcv -e -X$x1/$x2/$dx -Y$y1/$y2/$dy -Ostatic_greenfn.nc

grt static syn -Gstatic_greenfn.nc -Su1e12 -M$strike/$dip/$rake -N -e -Ostatic_syn.nc -s

grt okada -I6/3.464/2.7 -Su1e12 -Ds$depsrc -Dr$deprcv -M$strike/$dip/$rake -N -e -X$x1/$x2/$dx -Y$y1/$y2/$dy -Ookada.nc -s

cat > finite_faults.inp <<EOF
  #   X-start    Y-start     X-fin     Y-fin    Kode  shear(m)  reverse(m)  dip angle   top(km)   bot(km)
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000    -2.0000     0.0000     2.0000 100     0.0300      0.0400  60.0          2        4
EOF

grt static greenfn -Mhalfspace -Ds2/4/0.5 -Dr0 -e -R0/10/0.2 -Ofinite_greenfn.nc

grt static syn -Gfinite_greenfn.nc -Cfinite_faults.inp+i0.5/0.5 -N -e -X$x1/$x2/$dx -Y$y1/$y2/$dy -Ofinite_static_syn.nc -s

grt okada -I6/3.464/2.7 -Cfinite_faults.inp -Dr0 -N -e -X$x1/$x2/$dx -Y$y1/$y2/$dy -Ofinite_okada.nc -s

python plot.py

cp compare_displacement.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf *.nc halfspace finite_faults.inp
