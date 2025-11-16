#!/bin/bash

set -euo pipefail

grt lamb1 -h

grt lamb1 -P0.25 -T0/2/1e-3 -A30 > lamb1

python -u test_lamb1.py 

rm lamb1