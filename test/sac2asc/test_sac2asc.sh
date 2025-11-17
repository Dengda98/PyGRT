#!/bin/bash

set -euo pipefail

grt sac2asc -h

grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -OGRN
grt sac2asc GRN/milrow_2_3_10/VFZ.sac

rm -rf GRN