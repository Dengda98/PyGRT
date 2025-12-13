#!/bin/bash

set -euo pipefail

grt static greenfn -h
grt static_greenfn -h

grt static greenfn -M../milrow -D2/0 -X-3/3/0.2 -Y-2/2/0.2 -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -X-3/3/0.2 -Y-2/2/0.2 -e -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -X-3/3/0.2 -Y-2/2/0.2 -L20 -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -X-3/3/0.2 -Y-2/2/0.2 -K+k4+e1e-3 -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -X-3/3/0.2 -Y-2/2/0.2 -S -Ostgrn.nc

python -u test_static_greenfn.py

rm -rf stgrt* *.nc