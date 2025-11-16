#!/bin/bash

set -euo pipefail

grt static strain -h
grt static stress -h
grt static rotation -h


grt static greenfn -M../milrow -D2/0 -X-3/3/0.2 -Y-2/2/0.2 -e -Ostgrn.nc
grt static syn -S1e20 -Gstgrn.nc -e -Ostsyn.nc
grt static strain stsyn.nc
grt static stress stsyn.nc
grt static rotation stsyn.nc

grt static syn -S1e20 -Gstgrn.nc -e -N -Ostsyn_ZNE.nc
grt static strain stsyn_ZNE.nc
grt static stress stsyn_ZNE.nc
grt static rotation stsyn_ZNE.nc

python -u test_static_tensors.py

rm *.nc