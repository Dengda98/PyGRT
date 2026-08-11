#!/bin/bash

set -euo pipefail

grt eigenv -h

grt eigenv -M../milrow -SR -F0/1/0.01 -N -Cphase_R.nc
grt eigenv -M../milrow -SL -F0/1/0.01 -N -Cphase_L.nc
grt disp2asc -Cphase_R.nc -N > phase_R.txt

grt eigenv -M../milrow -SR -X1 > secfunc_R.txt
grt eigenv -M../milrow -SL -X1 > secfunc_L.txt

python -u test_eigenv.py

rm *.nc *.txt
