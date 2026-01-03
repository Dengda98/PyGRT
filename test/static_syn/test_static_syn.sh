#!/bin/bash

set -euo pipefail

grt static syn -h
grt static_syn -h

grt static greenfn -M../milrow -D2/0 -X-3/3/0.2 -Y-2/2/0.2 -e -Ostgrn.nc

grt static syn -S1e20 -Gstgrn.nc -Ostsyn.nc
grt static syn -S1e20 -F2/-1/4    -Gstgrn.nc -Ostsyn.nc
grt static syn -S1e20 -M77/88/111 -Gstgrn.nc -Ostsyn.nc
grt static syn -Su1e6 -M77/88/111 -Gstgrn.nc -Ostsyn.nc
grt static syn -Su1e6 -M77/88 -Gstgrn.nc -Ostsyn.nc
grt static syn -S1e20 -T1/-2/-5/0.5/3/1.2 -Gstgrn.nc -Ostsyn.nc

grt static syn -S1e20 -F2/-1/4 -e -Gstgrn.nc -Ostsyn.nc
grt static syn -S1e20 -F2/-1/4 -N -e -Gstgrn.nc -Ostsyn.nc

python -u test_static_syn.py

rm *.nc