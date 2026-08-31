#!/bin/bash

set -euo pipefail

grt lamb2 -h

grt lamb2 -P0.25 -T0/2/1e-3 -A30 -D60 > lamb2

python -u test_lamb2.py

rm -rf lamb2
