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

expect_fail "lamb3 source on the free surface" grt lamb3 -P0.25 -T0/0/1 -R10 -D0/1 -A0
expect_fail "lamb3 receiver on the free surface" grt lamb3 -P0.25 -T0/0/1 -R10 -D1/0 -A0
expect_fail "lamb3 zero horizontal distance" grt lamb3 -P0.25 -T0/0/1 -R0 -D2/1 -A0
expect_warn "lamb3 low Poisson ratio" "calculation is very likely to fail" grt lamb3 -P0.0005 -T0/0/1 -R10 -D2/1 -A0
expect_warn "lamb3 high Poisson ratio" "calculation is very likely to fail" grt lamb3 -P0.4995 -T0/0/1 -R10 -D2/1 -A0
expect_warn "lamb3 shallow source" "calculation is very likely to fail" grt lamb3 -P0.25 -T0/0/1 -R10 -D0.005/1 -A0
expect_warn "lamb3 small horizontal distance" "horizontal distance ratio" grt lamb3 -P0.25 -T0/0/1 -R1e-4 -D2/1 -A0

grt lamb3 -h
grt lamb3 -P0.25 -T0/2/1e-2 -R10 -D2/1 -S+slamb3_source+rlamb3_receiver -A30 > lamb3

python -u test_lamb3.py

rm -f lamb3 lamb3_source lamb3_receiver
