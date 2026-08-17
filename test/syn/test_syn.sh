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

grt syn -h

grt greenfn -M../milrow -D2/3 -N40/0.02 -R10 -OGRN_SINGLE -s
grt syn -GGRN_SINGLE -A22 -S1e20 -Osyn_single
grt syn -GGRN_SINGLE -Ds2 -Dr3 -R10 -A22 -S1e20 -Osyn_single_explicit
expect_fail "wrong single source depth is rejected" \
    grt syn -GGRN_SINGLE -Ds2.1 -A22 -S1e20 -Osyn_bad
expect_fail "wrong single receiver depth is rejected" \
    grt syn -GGRN_SINGLE -Dr3.1 -A22 -S1e20 -Osyn_bad
expect_fail "wrong single epicentral distance is rejected" \
    grt syn -GGRN_SINGLE -R11 -A22 -S1e20 -Osyn_bad

grt greenfn -M../milrow -Ds1,2 -Dr3 -N40/0.02 -R10 -OGRN_SRC_MULTI -s
grt syn -GGRN_SRC_MULTI -Ds1 -A22 -S1e20 -Osyn_src_multi
grt syn -GGRN_SRC_MULTI -Ds1 -Dr3 -A22 -S1e20 -Osyn_src_multi_explicit
expect_fail "multiple source depths require -Ds" \
    grt syn -GGRN_SRC_MULTI -A22 -S1e20 -Osyn_bad
expect_fail "wrong single receiver depth is rejected in source-depth library" \
    grt syn -GGRN_SRC_MULTI -Ds1 -Dr3.1 -A22 -S1e20 -Osyn_bad

grt greenfn -M../milrow -Ds2 -Dr1,3 -N40/0.02 -R10 -OGRN_RCV_MULTI -s
grt syn -GGRN_RCV_MULTI -Dr1 -A22 -S1e20 -Osyn_rcv_multi
grt syn -GGRN_RCV_MULTI -Ds2 -Dr1 -A22 -S1e20 -Osyn_rcv_multi_explicit
expect_fail "multiple receiver depths require -Dr" \
    grt syn -GGRN_RCV_MULTI -A22 -S1e20 -Osyn_bad
expect_fail "wrong single source depth is rejected in receiver-depth library" \
    grt syn -GGRN_RCV_MULTI -Dr1 -Ds2.1 -A22 -S1e20 -Osyn_bad

grt greenfn -M../milrow -D2/3 -N40/0.02 -R5,10 -OGRN_DIST_MULTI -s
grt syn -GGRN_DIST_MULTI -R10 -A22 -S1e20 -Osyn_dist_multi
expect_fail "multiple epicentral distances require -R" \
    grt syn -GGRN_DIST_MULTI -A22 -S1e20 -Osyn_bad

grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -e -OGRN
# 在同一根目录下准备多个深度组合和多个震中距，测试 syn 的精确检索
grt greenfn -M../milrow -D2/3 -N80/0.02 -R5 -OGRN -s
grt greenfn -M../milrow -D1/0 -N80/0.02 -R10 -OGRN -s

grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Osyn 
grt syn -GGRN -Ds1 -Dr0 -R10 -A22 -S1e20 -Osyn_multi_1_0_10
grt syn -GGRN -Ds2 -Dr3 -R10 -A22 -S1e20 -Osyn_multi_2_3_10
expect_fail "selectors cannot be used with a subdirectory" \
    grt syn -GGRN/milrow_2_3_10 -Ds2 -Dr3 -R10 -A22 -S1e20 -Osyn_bad
expect_fail "root directory requires exact matching distance" \
    grt syn -GGRN -Ds2 -Dr3 -R9 -A22 -S1e20 -Osyn_bad
expect_fail "root directory requires -R when depths match multiple distances" \
    grt syn -GGRN -Ds2 -Dr3 -A22 -S1e20 -Osyn_bad
grt syn -GGRN/milrow_2_3_10 -A22 -S1e16 -F-1/2/-4 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -M33/44/55 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -Su1e10 -M33/44/55 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -Su1e10 -M33/44 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20  -T1/-2/-5/0.5/3/1.2 -Osyn 

grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Dp/0.6 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Dt/0.2/0.4/0.7 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Dt/0.4/0.4/0.8 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Dr/1.2 -Osyn 
cat > tfile <<EOF
0       0.0
0.02    0.1
0.04    0.2
0.06    0.4
0.08    0.4
0.10    0.4
0.12    0.2
0.14    0.1
0.16    0.0
EOF
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -D0/tfile -Osyn 
rm -rf tfile

grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -M33/44/55 -I1 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -M33/44/55 -J1 -Osyn 

grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -N -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -e -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -N -e -Osyn 


python -u test_syn.py


rm -rf GRN GRN_SINGLE GRN_SRC_MULTI GRN_RCV_MULTI GRN_DIST_MULTI \
    syn syn_single syn_bad syn_src_multi syn_rcv_multi syn_dist_multi \
    syn_single_explicit syn_src_multi_explicit syn_rcv_multi_explicit \
    syn_subdir_bad syn_multi_1_0_10 syn_multi_2_3_10
