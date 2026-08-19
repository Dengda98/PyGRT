#!/bin/bash

set -euo pipefail

expect_fail() {
    # 期望命令失败；成功则报错退出
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
    # 期望命令成功且 stderr 含指定关键字
    local desc="$1"
    local key="$2"
    shift 2
    local tmp
    tmp=$(mktemp)
    set +e
    "$@" >"$tmp" 2>&1
    local ret=$?
    set -e
    if [ "$ret" -ne 0 ]; then
        echo "ERROR: expected success with warning but failed: $desc" >&2
        cat "$tmp" >&2
        rm -f "$tmp"
        exit 1
    fi
    if ! grep -q "$key" "$tmp"; then
        echo "ERROR: expected warning containing '$key': $desc" >&2
        cat "$tmp" >&2
        rm -f "$tmp"
        exit 1
    fi
    rm -f "$tmp"
    echo "OK (warned as expected): $desc"
}

grt static greenfn -h
grt static_greenfn -h

# -------------------- 单深度（旧 -D）--------------------
grt static greenfn -M../milrow -D2/0 -R0/4/0.2 -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/0.2 -e -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/0.2 -L20 -Ostgrn.nc

grt static greenfn -M../milrow -D0.2/0 -R0/4/0.2 -L20 -Cd -Ostgrn.nc
grt static greenfn -M../milrow -D0.2/0 -R0/4/0.2 -L20 -Cp -Ostgrn.nc
grt static greenfn -M../milrow -D0.2/0 -R0/4/0.2 -L20 -Cn -Ostgrn.nc

grt static greenfn -M../milrow -D2/0 -R0/4/0.2 -K+k4+e1e-3 -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/0.2 -S -Ostgrn.nc

# -X/-Y 仅保留一个二维网格输入测试
grt static greenfn -M../milrow -D2/0 -X-3/3/0.2 -Y-2/2/0.2 -Ostgrn_xy.nc

# -R
grt static greenfn -M../milrow -D2/0 -R0/10/0.1 -Ostgrn.nc
seq 0 0.1 10 > dists
grt static greenfn -M../milrow -D2/0 -Rdists -Ostgrn.nc
rm -rf dists
grt static greenfn -M../milrow -D2/0 -R2,3,5,8 -Ostgrn.nc

# boundary
grt static greenfn -M../milrow -D2/0 -R0/4/0.2 -BrF -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/0.2 -BhR -Ostgrn.nc
grt static greenfn -M../milrow -D2/0 -R0/4/0.2 -BrH -Ostgrn.nc

# -------------------- 多深度（-Ds / -Dr）--------------------
# 逗号列表
grt static greenfn -M../milrow -Ds1,2,3 -Dr0 -R0/4/1 -Ostgrn_multi.nc
# 等间距
grt static greenfn -M../milrow -Ds1/3/1 -Dr0/1/1 -R0/4/1 -Ostgrn_multi.nc
# 单深度也可用 -Ds/-Dr
grt static greenfn -M../milrow -Ds2 -Dr0 -R0/4/1 -Ostgrn_multi.nc
# 带位移偏导
grt static greenfn -M../milrow -Ds1,2 -Dr0,0.5 -R0/4/1 -e -Ostgrn_multi.nc
# 从文件读深度
printf "1\n2\n3\n" > depsrc_list
printf "0\n0.5\n" > deprcv_list
grt static greenfn -M../milrow -Dsdepsrc_list -Drdeprcv_list -R0/4/1 -Ostgrn_multi.nc
rm -f depsrc_list deprcv_list

# 多深度时 -S 应警告并忽略，但仍成功
expect_warn "multi-depth -S ignored" "ignored" \
    grt static greenfn -M../milrow -Ds1,2 -Dr0 -R0/4/1 -S -Ostgrn_multi.nc

# -------------------- 错误参数（CLI）--------------------
expect_fail "-D and -Ds/-Dr mutually exclusive" \
    grt static greenfn -M../milrow -D2/0 -Ds1,2 -Dr0 -R0/4/1 -Ostgrn_bad.nc

expect_fail "-Ds without -Dr" \
    grt static greenfn -M../milrow -Ds1,2 -R0/4/1 -Ostgrn_bad.nc

expect_fail "-Dr without -Ds" \
    grt static greenfn -M../milrow -Dr0 -R0/4/1 -Ostgrn_bad.nc

expect_fail "missing depth option" \
    grt static greenfn -M../milrow -R0/4/1 -Ostgrn_bad.nc

expect_fail "negative depth in -Ds" \
    grt static greenfn -M../milrow -Ds-1,2 -Dr0 -R0/4/1 -Ostgrn_bad.nc

expect_fail "negative depth in -D" \
    grt static greenfn -M../milrow -D-1/0 -R0/4/1 -Ostgrn_bad.nc

expect_fail "nonpositive depth spacing" \
    grt static greenfn -M../milrow -Ds1/3/0 -Dr0 -R0/4/1 -Ostgrn_bad.nc

expect_fail "depth start > end" \
    grt static greenfn -M../milrow -Ds3/1/1 -Dr0 -R0/4/1 -Ostgrn_bad.nc

expect_fail "non-ascending -R list" \
    grt static greenfn -M../milrow -D2/0 -R3,1,2 -Ostgrn_bad.nc

# -------------------- Python --------------------
python -u test_static_greenfn.py

rm -rf stgrt* *.nc
