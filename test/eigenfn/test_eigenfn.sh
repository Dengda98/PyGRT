#!/bin/bash

set -euo pipefail

grt eigenfn -h

grt eigenv -M../milrow -SR -F0/1/0.01 -N -Cphase_R.nc
grt eigenfn -Cphase_R.nc -F1 -N0/10 -Wegn_R.nc+z0/40/0.1 -Ugroup_R.nc -K+ccsens.nc+uusens.nc+z0.1+xegyint.nc
grt disp2asc -Ugroup_R.nc -N > group_R.txt


rm *.nc *.txt