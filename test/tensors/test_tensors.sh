#!/bin/bash

set -euo pipefail

grt strain -h
grt stress -h
grt rotation -h

grt greenfn -M../milrow -D2/0 -N600/0.02 -R10 -e -OGRN
grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -e -Osyn 
grt strain syn
grt stress syn
grt rotation syn

grt syn -GGRN/milrow_2_0_10 -A22 -S1e20 -e -N -Osyn_ZNE 
grt strain syn_ZNE 
grt stress syn_ZNE 
grt rotation syn_ZNE 

python -u test_tensors.py

rm GRN syn* -rf
