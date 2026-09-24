#!/bin/bash

set -euo pipefail

grt lamb1 -h

grt lamb1 -P0.25 -T0/2/1e-3 -A30 > lamb1
if grt lamb1 -P0.25 -T0/2/1e-3 -A30 -C0 > lamb1_invalid 2>&1; then
    echo "grt lamb1 should reject -C0." >&2
    exit 1
fi
grt lamb1 -P0.25 -T0/2/1e-3 -A30 -C2e-4 > lamb1_moving

python -u test_lamb1.py

rm -f lamb1 lamb1_invalid lamb1_moving
