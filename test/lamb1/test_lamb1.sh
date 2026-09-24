#!/bin/bash

set -euo pipefail

grt lamb1 -h

grt lamb1 -P0.25 -T0/2/1e-3 -A30 > lamb1
grt lamb1 -P0.25 -T0/2/1e-3 -A30 -QP,S > lamb1_PS
grt lamb1 -P0.25 -T0/2/1e-3 -A30 -QR > lamb1_R
if grt lamb1 -P0.25 -T0/2/1e-3 -A30 -C0 > lamb1_invalid 2>&1; then
    echo "grt lamb1 should reject -C0." >&2
    exit 1
fi
grt lamb1 -P0.25 -T0/2/1e-3 -A30 -C2e-4 > lamb1_moving

expect_warn() {
    local desc="$1"
    local key="$2"
    shift 2
    local output
    output=$("$@" 2>&1)
    if [[ "$output" != *"$key"* ]]; then
        echo "ERROR: expected warning containing '$key': $desc" >&2
        echo "$output" >&2
        exit 1
    fi
    echo "OK (warned as expected): $desc"
}

expect_warn "lamb1 invalid phase" "Available phases: P, S, R" grt lamb1 -P0.25 -T0/0/1 -A0 -QBAD,P
expect_warn "lamb1 slash phase separator" "Unsupported Lamb phase P/S" grt lamb1 -P0.25 -T0/0/1 -A0 -QP/S
expect_warn "lamb1 duplicated phase" "recorded only once" grt lamb1 -P0.25 -T0/0/1 -A0 -QP,P
expect_warn "lamb1 no valid phase" "output is all zeros" grt lamb1 -P0.25 -T0/0/1 -A0 -QBAD

python -u test_lamb1.py

rm -f lamb1 lamb1_invalid lamb1_moving lamb1_PS lamb1_R
