#!/bin/bash

set -euo pipefail

rm -rf GRN* *.nc

# BEGIN

M="milrow"
D="2/0"
R="100"
n="5"

# Calculate Full-wave
grt greenfn -M$M -D$D -R$R -N200/0.5+n$n -OGRN -e

# Calculate Surface-wave
grt eigenv -M$M -F0/1/0.01 -SR -N -Cphase_R.nc 
grt eigenv -M$M -F0/1/0.01 -SL -N -Cphase_L.nc 

# 0th order
grt modsum -Cphase_R.nc -D$D -R$R -N0 -OGRN_NM_0 -W$n -e
grt modsum -Cphase_L.nc -D$D -R$R -N0 -OGRN_NM_0 -W$n -e
# 1st order
grt modsum -Cphase_R.nc -D$D -R$R -N1 -OGRN_NM_1 -W$n -e
grt modsum -Cphase_L.nc -D$D -R$R -N1 -OGRN_NM_1 -W$n -e
# 2nd order
grt modsum -Cphase_R.nc -D$D -R$R -N2 -OGRN_NM_2 -W$n -e
grt modsum -Cphase_L.nc -D$D -R$R -N2 -OGRN_NM_2 -W$n -e
# all
grt modsum -Cphase_R.nc -D$D -R$R -N -OGRN_NM_all -W$n -e
grt modsum -Cphase_L.nc -D$D -R$R -N -OGRN_NM_all -W$n -e
# END

python compare_wave.py 0
python compare_wave.py 1
python compare_wave.py 2
python compare_wave.py all

rm -rf GRN* *.nc