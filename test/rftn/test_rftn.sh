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

grt rftn -h

# P 波水平射线参数，输出接收函数和绝对位移分量
grt rftn -M../milrow -P0.12 -TP -N64/0.05 -W -OC_P -s
test -f C_P/P_rftn.sac
test -f C_P/P_Z.sac
test -f C_P/P_R.sac

# SV 波入射角，指定反向入射层并覆盖可选频域参数
grt rftn -M../milrow -I20/2 -TS -N64/0.05+w0.9+n2+a+f -A0.8 -E0.5 -W -OC_S -s
test -f C_S/S_rftn.sac
test -f C_S/S_Z.sac
test -f C_S/S_R.sac

expect_fail "-P and -I are mutually exclusive" \
    grt rftn -M../milrow -P0.12 -I20 -TP -N8/0.1 -OC_BAD

expect_fail "one of -P and -I is required" \
    grt rftn -M../milrow -TP -N8/0.1 -OC_BAD

expect_fail "incidence angle must be below 90 degrees" \
    grt rftn -M../milrow -I90 -TP -N8/0.1 -OC_BAD

python -u test_rftn.py

rm -rf C_P C_S PY_P PY_S
