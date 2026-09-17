#!/bin/bash

set -euo pipefail

OUTPUT_DIRS=(lamb_surface lamb2_source lamb2_receiver lamb3_force lamb3_mt lamb_E_zero lamb_E_shift)
SURFACE_LOG=$(mktemp)
UNSUPPORTED_LOG=$(mktemp)
ZERO_ASC=$(mktemp)
SHIFT_ASC=$(mktemp)
trap 'rm -rf "${OUTPUT_DIRS[@]}" "$SURFACE_LOG" "$UNSUPPORTED_LOG" "$ZERO_ASC" "$SHIFT_ASC"' EXIT

grt lamb -h

# 地表源和地表台站只允许计算由 -F 指定的单力源
if grt lamb -H8.0/4.62/3.3 -N16/0.01 -R10 -Ds0 -Dr0 -A30 -S1 -O"${OUTPUT_DIRS[0]}" -s >"$UNSUPPORTED_LOG" 2>&1; then
    echo "surface source without -F unexpectedly succeeded" >&2
    exit 1
fi
grep -q "only the single force source specified by -F is supported" "$UNSUPPORTED_LOG"

grt lamb -H8.0/4.62/3.3 -N16/0.01 -R10 -Ds0 -Dr0 -A30 -M100/30/70 -S1e20 -O"${OUTPUT_DIRS[0]}" -s >"$UNSUPPORTED_LOG" 2>&1 && {
    echo "surface moment source unexpectedly succeeded" >&2
    exit 1
}
grep -q "only the single force source specified by -F is supported" "$UNSUPPORTED_LOG"

grt lamb -H8.0/4.62/3.3 -N16/0.01 -R10 -Ds0 -Dr0 -A30 -F0.5/-1/2 -S1 -O"${OUTPUT_DIRS[0]}" -e -n -s >"$SURFACE_LOG" 2>&1
grep -q -- "-e is ignored" "$SURFACE_LOG"
test -f "${OUTPUT_DIRS[0]}/Z.sac"
test -f "${OUTPUT_DIRS[0]}/N.sac"
test -f "${OUTPUT_DIRS[0]}/E.sac"
test ! -e "${OUTPUT_DIRS[0]}/zZ.sac"
test ! -e "${OUTPUT_DIRS[0]}/nZ.sac"
test ! -e "${OUTPUT_DIRS[0]}/eZ.sac"

# 第二类 Lamb：地下源、地表台站，以及地表源、地下台站
grt lamb -H8.0/4.62/3.3 -N16/0.01 -R10 -Ds5 -Dr0 -A30 -M100/30/70 -S1e20 -O"${OUTPUT_DIRS[1]}" -e -n -s
grt lamb -H8.0/4.62/3.3 -N16/0.01 -R10 -Ds0 -Dr5 -A30 -T1/0/0/1/0/1 -S1e20 -O"${OUTPUT_DIRS[2]}" -s
test -f "${OUTPUT_DIRS[1]}/Z.sac"
test -f "${OUTPUT_DIRS[1]}/N.sac"
test -f "${OUTPUT_DIRS[1]}/E.sac"
test -f "${OUTPUT_DIRS[1]}/zZ.sac"
test -f "${OUTPUT_DIRS[1]}/nZ.sac"
test -f "${OUTPUT_DIRS[1]}/eZ.sac"
test -f "${OUTPUT_DIRS[2]}/Z.sac"

# 第三类 Lamb：单力源、矩张量源、时间函数、积分和微分
grt lamb -H8.0/4.62/3.3 -N16/0.01 -R10 -Ds5 -Dr1 -A30 -F0.5/-1/2 -S1e15 -O"${OUTPUT_DIRS[3]}" -e -s
grt lamb -H8.0/4.62/3.3 -N16/0.01 -R10 -Ds5 -Dr1 -A30 -T1/0/0/1/0/1 -S1e20 -Dp/0.05 -I1 -J1 -O"${OUTPUT_DIRS[4]}" -n -s
test -f "${OUTPUT_DIRS[3]}/Z.sac"
test -f "${OUTPUT_DIRS[3]}/R.sac"
test -f "${OUTPUT_DIRS[3]}/zZ.sac"
test -f "${OUTPUT_DIRS[3]}/rZ.sac"
test -f "${OUTPUT_DIRS[3]}/tZ.sac"
test -f "${OUTPUT_DIRS[4]}/Z.sac"
test -f "${OUTPUT_DIRS[4]}/E.sac"
test -f "${OUTPUT_DIRS[4]}/sig.sac"

# -E 必须同时改变求解器的实际时间序列和 SAC 起始时刻
grt lamb -H8.0/4.62/3.3 -N500/0.01 -R10 -Ds5 -Dr1 -A30 -F0.5/-1/2 -S1e20 -O"${OUTPUT_DIRS[5]}" -s
grt lamb -H8.0/4.62/3.3 -N500/0.01 -R10 -Ds5 -Dr1 -A30 -F0.5/-1/2 -S1e20 -E0.2 -O"${OUTPUT_DIRS[6]}" -s
grt sac2asc "${OUTPUT_DIRS[5]}/Z.sac" >"$ZERO_ASC"
grt sac2asc "${OUTPUT_DIRS[6]}/Z.sac" >"$SHIFT_ASC"
awk 'NR == 1 && ($1 < 0.199 || $1 > 0.201) { exit 1 }' "$SHIFT_ASC"
awk '
    NR == FNR {
        zero[FNR] = $2
        if ($2 > 1e-6 || $2 < -1e-6) {
            has_signal = 1
        }
        next
    }
    FNR <= 479 {
        difference = $2 - zero[FNR + 20]
        if (difference > 1e-5 || difference < -1e-5) {
            exit 1
        }
    }
    END {
        if (!has_signal) {
            exit 1
        }
    }
' "$ZERO_ASC" "$SHIFT_ASC"

python -u test_lamb.py
python -u compare_lamb_syn.py
