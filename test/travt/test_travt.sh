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

grt travt -h 

grt travt -M../milrow -D2/0 -R2,3,4,5
cat > dists <<EOF
2
3
4
5
EOF
grt travt -M../milrow -D2/0 -Rdists
rm -rf dists

expect_fail "non-ascending -R list" \
    grt travt -M../milrow -D2/0 -R3,1,2

python -u test_travt.py
