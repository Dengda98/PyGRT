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

grt static syn -h
grt static_syn -h

# -------------------- 单深度 XY 库 --------------------
grt static greenfn -M../milrow -D2/0 -X-3/3/0.2 -Y-2/2/0.2 -e -Ostgrn.nc

grt static syn -S1e20 -Gstgrn.nc -Ostsyn.nc
grt static syn -S1e20 -F2/-1/4    -Gstgrn.nc -Ostsyn.nc
grt static syn -S1e20 -M77/88/111 -Gstgrn.nc -Ostsyn.nc
grt static syn -Su1e6 -M77/88/111 -Gstgrn.nc -Ostsyn.nc
grt static syn -Su1e6 -M77/88 -Gstgrn.nc -Ostsyn.nc
grt static syn -S1e20 -T1/-2/-5/0.5/3/1.2 -Gstgrn.nc -Ostsyn.nc

grt static syn -S1e20 -F2/-1/4 -e -Gstgrn.nc -Ostsyn.nc
grt static syn -S1e20 -F2/-1/4 -N -e -Gstgrn.nc -Ostsyn.nc

# 合成时指定新 XY 网格（震中距插值）
grt static syn -S1e20 -Gstgrn.nc -X-2/2/0.5 -Y-1/1/0.5 -Ostsyn_xy.nc

# -------------------- -R 建库后合成 --------------------
grt static greenfn -M../milrow -D2/0 -R0/8/1 -e -Ostgrn_r.nc
# 延用库的震中距轴
grt static syn -S1e20 -Gstgrn_r.nc -Ostsyn_r.nc
# 指定新 XY 网格
grt static syn -S1e20 -Gstgrn_r.nc -X-3/3/1 -Y-2/2/1 -e -Ostsyn_rxy.nc

# -------------------- 多深度库：点源 -Ds；深度插值；-Q --------------------
grt static greenfn -M../milrow -Ds1,2,3 -Dr0 -R0/8/1 -e -Ostgrn_md.nc
grt static syn -Gstgrn_md.nc -Su1e16 -Ds2 -Ostsyn_md.nc
# 震源深度插值（库节点之间）
grt static syn -Gstgrn_md.nc -Su1e16 -Ds1.5 -X-2/2/1 -Y-2/2/1 -e -Ostsyn_interp.nc

cat > rcv_pts.txt <<'EOF'
# north east depth (km)
0 0 0
1 2 0
-1 1 0
EOF
grt static syn -Gstgrn_md.nc -Su1e16 -Ds2 -Qrcv_pts.txt -Ostsyn_q.nc

# 多台站深度：必须 -Dr
grt static greenfn -M../milrow -Ds2 -Dr0,0.5 -R0,5 -e -Ostgrn_mr.nc
grt static syn -Gstgrn_mr.nc -S1e20 -Dr0.25 -X-2/2/1 -Y-2/2/1 -Ostsyn_dr.nc

# -------------------- 有限断层（库震源深度覆盖断层 top/bot）--------------------
# W=(2.8-1.2)/sin(90°)=1.6 km，dW=1 → 末块短于 dW，用于覆盖余数子断层中心
cat > cfaults_tiny.inp <<'EOF'
  #   X-start    Y-start     X-fin     Y-fin    Kode  shear(m)  reverse(m)  dip angle   top(km)   bot(km)
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000     0.0000     2.0000     0.0000 100 0.1000     0.0000     90.00         1.2000    2.8000
EOF
grt static syn -Gstgrn_md.nc -Ccfaults_tiny.inp+i1/1 -Ostsyn_ff.nc
grt static syn -Gstgrn_md.nc -Ccfaults_tiny.inp+i1/1 -e -Qrcv_pts.txt -Ostsyn_ffq.nc

# -------------------- 错误参数 --------------------
expect_fail "multi-src library requires -Ds" \
    grt static syn -Gstgrn_md.nc -S1e20 -Ostsyn_bad.nc

expect_fail "multi-rcv library requires -Dr" \
    grt static syn -Gstgrn_mr.nc -S1e20 -Ostsyn_bad.nc

expect_fail "single-rcv library forbids -Dr" \
    grt static syn -Gstgrn.nc -S1e20 -Dr0 -Ostsyn_bad.nc

expect_fail "-Q mutually exclusive with -X/-Y" \
    grt static syn -Gstgrn_md.nc -S1e20 -Ds2 -Qrcv_pts.txt -X-1/1/1 -Y-1/1/1 -Ostsyn_bad.nc

expect_fail "-Q mutually exclusive with -Dr" \
    grt static syn -Gstgrn_md.nc -S1e20 -Ds2 -Dr0 -Qrcv_pts.txt -Ostsyn_bad.nc

expect_fail "finite fault requires ndepsrc>1" \
    grt static syn -Gstgrn.nc -Ccfaults_tiny.inp -Ostsyn_bad.nc

expect_fail "finite fault mutually exclusive with -S" \
    grt static syn -Gstgrn_md.nc -S1e20 -Ccfaults_tiny.inp -Ostsyn_bad.nc

cat > cfaults_bad_dip.inp <<'EOF'
  #   X-start    Y-start     X-fin     Y-fin    Kode  shear(m)  reverse(m)  dip angle   top(km)   bot(km)
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000     0.0000     2.0000     0.0000 100 0.1000     0.0000      0.00         1.2000    2.8000
EOF
expect_fail "finite fault dip must be in (0, 90]" \
    grt static syn -Gstgrn_md.nc -Ccfaults_bad_dip.inp -Ostsyn_bad.nc

cat > cfaults_bad_bot.inp <<'EOF'
  #   X-start    Y-start     X-fin     Y-fin    Kode  shear(m)  reverse(m)  dip angle   top(km)   bot(km)
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000     0.0000     2.0000     0.0000 100 0.1000     0.0000     90.00         2.8000    1.2000
EOF
expect_fail "finite fault bot must be greater than top" \
    grt static syn -Gstgrn_md.nc -Ccfaults_bad_bot.inp -Ostsyn_bad.nc

python -u test_static_syn.py

rm -rf *.nc rcv_pts.txt cfaults_tiny.inp cfaults_bad_dip.inp cfaults_bad_bot.inp
