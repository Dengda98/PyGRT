#!/bin/bash

set -euo pipefail

expect_fail() {
    local desc="$1"
    shift
    set +e
    "$@" >/dev/null 2>&1
    local ret=$?
    set -e
    if [ "$ret" -eq 0 ]; then
        echo "ERROR: expected failure but succeeded: $desc" >&2
        exit 1
    fi
    echo "OK (failed as expected): $desc"
}

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

expect_fail "lamb2 missing depth option" grt lamb2 -P0.25 -T0/0/1 -R10 -A0
expect_fail "lamb2 both -Ds and -Dr" grt lamb2 -P0.25 -T0/0/1 -R10 -Ds5 -Dr3 -A0
expect_fail "lamb2 old -D option" grt lamb2 -P0.25 -T0/0/1 -R10 -D5 -A0
expect_fail "lamb2 source on the free surface" grt lamb2 -P0.25 -T0/0/1 -R10 -Ds0 -A0
expect_fail "lamb2 receiver on the free surface" grt lamb2 -P0.25 -T0/0/1 -R10 -Dr0 -A0
expect_fail "lamb2 negative source depth" grt lamb2 -P0.25 -T0/0/1 -R10 -Ds-1 -A0
expect_fail "lamb2 negative receiver depth" grt lamb2 -P0.25 -T0/0/1 -R10 -Dr-1 -A0
expect_fail "lamb2 negative horizontal distance" grt lamb2 -P0.25 -T0/0/1 -R-1 -Ds5 -A0
expect_fail "lamb2 zero horizontal distance" grt lamb2 -P0.25 -T0/0/1 -R0 -Ds5 -A0
expect_warn "lamb2 low Poisson ratio" "calculation is very likely to fail" grt lamb2 -P0.0005 -T0/0/1 -R10 -Ds5 -A0
expect_warn "lamb2 high Poisson ratio" "calculation is very likely to fail" grt lamb2 -P0.4995 -T0/0/1 -R10 -Ds5 -A0
expect_warn "lamb2 shallow source" "calculation is very likely to fail" grt lamb2 -P0.25 -T0/0/1 -R10 -Ds0.005 -A0
expect_warn "lamb2 shallow receiver" "calculation is very likely to fail" grt lamb2 -P0.25 -T0/0/1 -R10 -Dr0.005 -A0
expect_warn "lamb2 small horizontal distance" "horizontal distance ratio" grt lamb2 -P0.25 -T0/0/1 -R1e-4 -Ds5 -A0

grt lamb2 -h

grt lamb2 -P0.25 -T0/2/1e-3 -R10 -Ds5 -S+slamb2_source+rlamb2_receiver -A30 > lamb2
grt lamb2 -P0.25 -T0/2/1e-3 -R10 -Dr5 -S+slamb2_surface_source+rlamb2_surface_receiver -A30 > lamb2_surface

python -u test_lamb2.py

rm -f lamb2 lamb2_source lamb2_receiver lamb2_surface lamb2_surface_source lamb2_surface_receiver
