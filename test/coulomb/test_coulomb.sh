#!/bin/bash

set -euo pipefail

expect_fail() {
    local desc="$1"
    local log="$2"
    shift 2
    if "$@" >"$log" 2>&1; then
        echo "ERROR: expected failure but succeeded: $desc" >&2
        exit 1
    fi
    echo "OK (failed as expected): $desc"
}

grt static coulomb -h
grt static_coulomb -h

cat > rcv_pts.txt <<'EOF'
# north east depth (km)
0 0 0
1 2 0
-1 1 0
EOF

cat > rcv_faults_defined.inp <<'EOF'
  #   X-start    Y-start      X-fin      Y-fin   Kode  value1      value2       dip       top       bot
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000     0.0000     2.0000     0.0000   100     0.1000      0.1000     90.00       1.0000     3.0000
EOF

# 使用同一个多接收深度格林函数，覆盖 grid、普通 points 和有限接收断层
grt static greenfn -M../milrow -Ds2 -Dr0/3/1 -R0/8/1 -e -Ostgrn.nc

# -------------------- grid：ZRT --------------------
grt static syn -Gstgrn.nc -S1e20 -Ds2 -Dr0 -e -X-2/2/1 -Y-1/1/1 -Ogrid_zrt.nc
grt static stress grid_zrt.nc
grt static sproj -Ggrid_zrt.nc -M33/44/55
ncdump -v sigma_n,tau_s grid_zrt.nc | sed -n '/^data:/,$p' > grid_zrt_before.txt
expect_fail "missing friction coefficient" grid_no_f.log \
    grt static coulomb -Ggrid_zrt.nc
expect_fail "negative friction coefficient" grid_negative_f.log \
    grt static coulomb -Ggrid_zrt.nc -F-0.1
expect_fail "invalid friction coefficient" grid_invalid_f.log \
    grt static coulomb -Ggrid_zrt.nc -Fabc
grt static coulomb -Ggrid_zrt.nc -F0.6
ncdump -v sigma_n,tau_s grid_zrt.nc | sed -n '/^data:/,$p' > grid_zrt_after.txt
cmp grid_zrt_before.txt grid_zrt_after.txt
ncdump -h grid_zrt.nc | rg 'coulomb\(north, east\)' >/dev/null

# 再次运行检查已有 coulomb 变量的覆盖警告
grt static coulomb -Ggrid_zrt.nc -F0.8 > grid_overwrite.log 2>&1
rg 'already exists and will be overwritten' grid_overwrite.log >/dev/null

# -------------------- grid：ZNE --------------------
grt static syn -Gstgrn.nc -S1e20 -Ds2 -Dr0 -N -e -X-2/2/1 -Y-1/1/1 -Ogrid_zne.nc
grt static stress grid_zne.nc
grt static sproj -Ggrid_zne.nc -M33/44/55
grt static coulomb -Ggrid_zne.nc -F0.6
ncdump -h grid_zne.nc | rg 'coulomb\(north, east\)' >/dev/null

# -------------------- 普通 points：ZRT --------------------
grt static syn -Gstgrn.nc -S1e20 -Ds2 -Qrcv_pts.txt -e -Opoints_zrt.nc
grt static stress points_zrt.nc
grt static sproj -Gpoints_zrt.nc -M33/44/55
grt static coulomb -Gpoints_zrt.nc -F0.25
ncdump -h points_zrt.nc | rg 'coulomb\(point\)' >/dev/null

# -------------------- 有限接收断层：ZRT --------------------
grt static syn -Gstgrn.nc -S1e20 -Ds2 -Rrcv_faults_defined.inp+i1/1 -e -Ofinite_zrt.nc
grt static stress finite_zrt.nc
grt static sproj -Gfinite_zrt.nc
grt static coulomb -Gfinite_zrt.nc -F0.4
ncdump -h finite_zrt.nc | rg 'coulomb\(point\)' >/dev/null

conda run -n pygrt-dev python -u test_coulomb.py

rm -f rcv_pts.txt rcv_faults_defined.inp stgrn.nc \
    grid_zrt.nc grid_zne.nc points_zrt.nc finite_zrt.nc \
    grid_zrt_before.txt grid_zrt_after.txt grid_no_f.log \
    grid_negative_f.log grid_invalid_f.log grid_overwrite.log

echo "test_coulomb.sh: all checks passed"
