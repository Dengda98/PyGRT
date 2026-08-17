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

# 计算面波频散, eigenv -F 和 greenfn -N 相匹配
grt eigenv -M../milrow -F0/1/0.01 -SR -N -Cphase_R.nc
grt eigenv -M../milrow -F0/1/0.01 -SL -N -Cphase_L.nc

# 模态叠加得到格林函数
# 仅 0 阶
grt modsum -Cphase_R.nc -D2/1 -R100 -N0 -OGRN_NM_0 -W5 -e
# 输出根目录应保留模型文件副本（路径来自频散文件中记录的模型）
test -f GRN_NM_0/milrow
grt modsum -Cphase_L.nc -D2/1 -R100 -N0 -OGRN_NM_0 -W5 -e
# 仅 1 阶
grt modsum -Cphase_R.nc -D2/1 -R100 -N1 -OGRN_NM_1 -W5 -e
grt modsum -Cphase_L.nc -D2/1 -R100 -N1 -OGRN_NM_1 -W5 -e
# 仅 2 阶
grt modsum -Cphase_R.nc -D2/1 -R100 -N2 -OGRN_NM_2 -W5 -e
grt modsum -Cphase_L.nc -D2/1 -R100 -N2 -OGRN_NM_2 -W5 -e
# 全部
grt modsum -Cphase_R.nc -D2/1 -R100 -N -OGRN_NM_all -W5 -e
grt modsum -Cphase_L.nc -D2/1 -R100 -N -OGRN_NM_all -W5 -e

# 多震源/台站深度
grt modsum -Cphase_R.nc -Ds1,2 -Dr0,1 -R100 -N0 -OGRN_NM_MULTI -W2 -e
test -f GRN_NM_MULTI/milrow_1_0_100/EXZ.sac
test -f GRN_NM_MULTI/milrow_2_1_100/EXZ.sac

expect_fail "-Ds without -Dr" \
    grt modsum -Cphase_R.nc -Ds1,2 -R100 -N0 -OGRN_bad

expect_fail "-D and -Ds/-Dr are mutually exclusive" \
    grt modsum -Cphase_R.nc -D2/1 -Ds1,2 -Dr0 -R100 -N0 -OGRN_bad

expect_fail "non-ascending -R list" \
    grt modsum -Cphase_R.nc -D2/1 -R100,50 -N0 -OGRN_bad

rm -rf GRN* *.nc
