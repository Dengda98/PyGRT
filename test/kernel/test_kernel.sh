#!/bin/bash

set -euo pipefail

grt kernel -h

grt kernel -M../thin1 -D0.03/0 -F0/25/1 -C0.01 -OKERNEL
grt kernel -M../seafloor -D5/0.1 -F0/1/1e-3 -C0.01 -OKERNEL

rm KERNEL -rf