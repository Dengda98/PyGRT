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

grt greenfn -h

grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -OGRN
# 输出根目录应保留模型文件副本
test -f GRN/milrow
grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -e -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02+w0.6+n10 -R10 -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02 -H1/10 -R10 -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02 -L20 -R10 -OGRN

grt greenfn -M../milrow -D2/3 -N600/0.02 -K+k4+s1.2+e1e-3+v1.5 -R10 -OGRN

grt greenfn -M../milrow -D2/0 -N1400/1 -L+a1e-3 -E-2/9 -R2000 -OGRN
grt greenfn -M../milrow -D2/0 -N1400/1+a -L+a1e-3 -E-2/9 -R2000 -OGRN
grt greenfn -M../milrow -D2/0 -N1400/1 -L+a1e-3 -Ep-20 -R2000 -OGRN

grt greenfn -M../milrow -D0.1/0 -N600/0.02 -R10 -OGRN -Cd
grt greenfn -M../milrow -D0.1/0 -N600/0.02 -R10 -OGRN -Cp
grt greenfn -M../milrow -D0.1/0 -N600/0.02 -R10 -OGRN -Cn

grt greenfn -M../milrow -D2/3 -N600/0.02 -L20 -R10 -S -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02 -L20 -R10 -S1,10,20 -OGRN

# boundary
grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -BrF -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -BhR -OGRN
grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -BrH -OGRN


# multi distance
grt greenfn -M../milrow -D2/0 -N600/0.02 -R6,8,10 -OGRN
cat > dists <<EOF
6
8
10
EOF
grt greenfn -M../milrow -D2/0 -N600/0.02 -Rdists -OGRN
rm -rf dists
grt greenfn -M../milrow -D2/0 -N600/0.02 -R6/10/2 -OGRN
test -f GRN/milrow_2_0_6/EXZ.sac
test -f GRN/milrow_2_0_8/EXZ.sac
test -f GRN/milrow_2_0_10/EXZ.sac

printf '6\n8\n10' > dists_no_newline
grt greenfn -M../milrow -D2/0 -N600/0.02 -Rdists_no_newline -OGRN
rm -f dists_no_newline

# multi source/receiver depths
grt greenfn -M../milrow -Ds1,2 -Dr0,1 -N80/0.02 -R5 -OGRN_MULTI -s
test -f GRN_MULTI/milrow_1_0_5/EXZ.sac
test -f GRN_MULTI/milrow_2_1_5/EXZ.sac

expect_fail "non-ascending -R list" \
    grt greenfn -M../milrow -D2/0 -N600/0.02 -R3,1,2 -OGRN_bad

expect_fail "-D and -Ds/-Dr are mutually exclusive" \
    grt greenfn -M../milrow -D2/0 -Ds1,2 -Dr0 -N80/0.02 -R5 -OGRN_bad

expect_fail "-Ds without -Dr" \
    grt greenfn -M../milrow -Ds1,2 -N80/0.02 -R5 -OGRN_bad


python -u test_greenfn.py

rm -rf GRN
rm -rf GRN_MULTI
rm -rf GRN_grtstats
