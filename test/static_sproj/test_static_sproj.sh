#!/bin/bash

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

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

grt static sproj -h
grt static_sproj -h

cat > rcv_pts.txt <<'EOF'
# north east depth (km)
0 0 0
1 2 0
-1 1 0
EOF

cat > rcv_pts_6.txt <<'EOF'
# north east depth strike dip rake
0 0 0 10 20 30
1 2 0 40 50 60
-1 1 0 70 80 90
EOF

cat > rcv_pts_new_6.txt <<'EOF'
# north east depth strike dip rake
0 0 0 120 35 -40
1 2 0 210 55 80
-1 1 0 300 25 140
EOF

cat > rcv_pts_reordered_6.txt <<'EOF'
# north east depth strike dip rake
1 2 0 120 35 -40
0 0 0 210 55 80
-1 1 0 300 25 140
EOF

cat > rcv_faults_undefined.inp <<'EOF'
  #   X-start    Y-start      X-fin      Y-fin   Kode  value1      value2       dip       top       bot
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000     0.0000     2.0000     0.0000   100     0.0000      0.0000     90.00       1.0000     3.0000
  2     0.0000     3.0000     2.0000     3.0000   200     0.1000      0.1000     60.00       1.0000     3.0000
EOF

cat > rcv_faults_defined.inp <<'EOF'
  #   X-start    Y-start      X-fin      Y-fin   Kode  value1      value2       dip       top       bot
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx  xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
  1     0.0000     0.0000     2.0000     0.0000   100     0.1000      0.1000     90.00       1.0000     3.0000
EOF

python - <<'PY'
import numpy as np
from scipy.io import netcdf_file

with netcdf_file("finite_no_rake.nc", mode="w") as nc:
    nc.createDimension("point", 2)
    nc.createDimension("nfault", 1)

    point_data = {
        "north": (0.0, 1.0),
        "east": (0.0, 1.0),
        "depth": (1.0, 1.0),
        "stress_ZZ": (1.0, 2.0),
        "stress_ZN": (0.0, 0.0),
        "stress_ZE": (0.0, 0.0),
        "stress_NN": (3.0, 4.0),
        "stress_NE": (0.0, 0.0),
        "stress_EE": (5.0, 6.0),
    }
    for name, values in point_data.items():
        variable = nc.createVariable(name, "d", ("point",))
        variable[:] = np.asarray(values, dtype=np.float64)

    for name, value in {"strike": 20.0, "dip": 45.0}.items():
        variable = nc.createVariable(name, "d", ("nfault",))
        variable[:] = value

    offset = nc.createVariable("offset", "i", ("nfault",))
    offset[:] = 2
    nc.layout = "points"
    nc.rot2ZNE = 1
PY

# 统一使用包含多个接收深度的库，覆盖普通点和有限接收断层
grt static greenfn -M../milrow -Ds2 -Dr0/3/1 -R0/8/1 -e -Ostgrn.nc

expect_fail "missing finite rake variable requires -M" finite_no_rake_no_m.log \
    grt static sproj -Gfinite_no_rake.nc
grt static sproj -Gfinite_no_rake.nc -M40

# -------------------- grid：手动三要素，ZRT/ZNE --------------------
grt static syn -Gstgrn.nc -S1e20 -Ds2 -Dr0 -e -X-2/2/1 -Y-1/1/1 -Ogrid_zrt.nc
grt static stress grid_zrt.nc
cp grid_zrt.nc grid_before.nc
expect_fail "grid requires -M" grid_no_m.log grt static sproj -Ggrid_zrt.nc
grt static sproj -Ggrid_zrt.nc -M33/44/55
python - <<'PY'
import numpy as np
from scipy.io import netcdf_file

with netcdf_file("grid_before.nc", mode="r", mmap=False) as before, netcdf_file(
    "grid_zrt.nc", mode="r", mmap=False
) as after:
    for name in ("north", "east", "stress_ZZ"):
        np.testing.assert_array_equal(before.variables[name][:], after.variables[name][:])
    assert after.variables["sigma_n"].dimensions == ("north", "east")
    assert after.variables["tau_s"].dimensions == ("north", "east")
PY

# 同一个文件再次运行，检查已有变量覆盖警告
grt static sproj -Ggrid_zrt.nc -M33/44/55 > grid_overwrite.log 2>&1
grep -Fq 'already exists and will be overwritten' grid_overwrite.log

grt static syn -Gstgrn.nc -S1e20 -Ds2 -Dr0 -N -e -X-2/2/1 -Y-1/1/1 -Ogrid_zne.nc
grt static stress grid_zne.nc
grt static sproj -Ggrid_zne.nc -M33/44/55

# -------------------- 普通 points：无几何、文件几何、-Q 新几何 --------------------
grt static syn -Gstgrn.nc -S1e20 -Ds2 -Qrcv_pts.txt -e -Opoints_plain.nc
grt static stress points_plain.nc
expect_fail "plain points require -M or -Q" points_plain_no_geometry.log \
    grt static sproj -Gpoints_plain.nc
grt static sproj -Gpoints_plain.nc -M33/44/55
cp points_plain.nc points_plain_manual.nc

grt static syn -Gstgrn.nc -S1e20 -Ds2 -Qrcv_pts_6.txt -e -Opoints_geometry.nc
grt static stress points_geometry.nc
grt static sproj -Gpoints_geometry.nc
cp points_geometry.nc points_geometry_from_file.nc
grt static sproj -Gpoints_geometry.nc -M120/35/-40 > points_manual.log 2>&1
grep -Fq 'manual geometry' points_manual.log

grt static sproj -Gpoints_plain.nc -Qrcv_pts_new_6.txt > points_Q.log 2>&1
grep -Fq 'already exists and will be overwritten' points_Q.log
expect_fail "-Q requires six columns" points_Q_3cols.log \
    grt static sproj -Gpoints_plain.nc -Qrcv_pts.txt
expect_fail "-Q point order must match" points_Q_order.log \
    grt static sproj -Gpoints_plain.nc -Qrcv_pts_reordered_6.txt

# -------------------- 有限接收断层：未定义与已定义 rake --------------------
grt static syn -Gstgrn.nc -S1e20 -Ds2 -Rrcv_faults_undefined.inp+i1/1 -e -Ofinite_undefined.nc
grt static stress finite_undefined.nc
python - <<'PY'
from scipy.io import netcdf_file

with netcdf_file("finite_undefined.nc", mode="r", mmap=False) as nc:
    assert "nfault" in nc.dimensions
PY
expect_fail "finite receiver with undefined rake requires -M" finite_no_m.log \
grt static sproj -Gfinite_undefined.nc
grt static sproj -Gfinite_undefined.nc -M55
cp finite_undefined.nc finite_undefined_partial.nc
grt static sproj -Gfinite_undefined.nc -M55+f > finite_force.log 2>&1
grep -Fq 'already exists and will be overwritten' finite_force.log

grt static syn -Gstgrn.nc -S1e20 -Ds2 -Rrcv_faults_defined.inp+i1/1 -e -Ofinite_defined.nc
grt static stress finite_defined.nc
grt static sproj -Gfinite_defined.nc
expect_fail "defined finite rake rejects non-forcing -M" finite_defined_m.log \
    grt static sproj -Gfinite_defined.nc -M55
grt static sproj -Gfinite_defined.nc -M55+f

python -u test_sproj.py

rm -f rcv_pts.txt rcv_pts_6.txt rcv_pts_new_6.txt rcv_pts_reordered_6.txt \
    rcv_faults_undefined.inp rcv_faults_defined.inp stgrn.nc \
    grid_zrt.nc grid_zne.nc points_plain.nc points_plain_manual.nc \
    points_geometry.nc points_geometry_from_file.nc finite_no_rake.nc finite_undefined.nc \
    finite_undefined_partial.nc finite_defined.nc \
    grid_before.nc grid_no_m.log grid_overwrite.log \
    points_plain_no_geometry.log points_manual.log points_Q.log points_Q_3cols.log \
    points_Q_order.log finite_no_rake_no_m.log finite_no_m.log finite_force.log finite_defined_m.log

echo "test_sproj.sh: all checks passed"
