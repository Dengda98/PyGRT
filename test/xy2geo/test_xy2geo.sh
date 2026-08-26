#!/bin/bash

set -euo pipefail

grt xy2geo -h
grt geo2xy -h

python - <<'PY'
import numpy as np
from scipy.io import netcdf_file


with netcdf_file("grid_local.nc", mode="w") as nc:
    nc.createDimension("north", 3)
    nc.createDimension("east", 3)
    north = nc.createVariable("north", "d", ("north",))
    east = nc.createVariable("east", "d", ("east",))
    field = nc.createVariable("field", "d", ("north", "east"))
    north[:] = [-10.0, 0.0, 10.0]
    east[:] = [-200.0, 0.0, 200.0]
    field[:] = np.arange(9, dtype=np.float64).reshape(3, 3)
    nc.layout = "grid"
    nc.computeType = "test"


with netcdf_file("points_local.nc", mode="w") as nc:
    nc.createDimension("point", 3)
    for name, values in {
        "north": [0.0, 10.0, -20.0],
        "east": [-200.0, 0.0, 200.0],
        "depth": [0.0, 1.0, 2.0],
        "field": [1.0, 2.0, 3.0],
    }.items():
        variable = nc.createVariable(name, "d", ("point",))
        variable[:] = values
    nc.layout = "points"
    nc.computeType = "test"


with netcdf_file("finite_local.nc", mode="w") as nc:
    nc.createDimension("point", 4)
    nc.createDimension("nfault", 2)
    for name, values in {
        "north": [0.0, 10.0, 20.0, 30.0],
        "east": [-2.0, -1.0, 1.0, 2.0],
        "depth": [1.0, 1.0, 2.0, 2.0],
        "field": [4.0, 5.0, 6.0, 7.0],
    }.items():
        variable = nc.createVariable(name, "d", ("point",))
        variable[:] = values
    for name, values in {
        "strike": [10.0, 20.0],
        "dip": [30.0, 40.0],
        "rake": [50.0, 60.0],
    }.items():
        variable = nc.createVariable(name, "d", ("nfault",))
        variable[:] = values
    offset = nc.createVariable("offset", "i", ("nfault",))
    offset[:] = [2, 4]
    nc.layout = "points"
    nc.computeType = "test"
PY

cat > local_points.txt <<'EOF'
# north east station

  # an indented comment
  0.0    0.0 station_A extra text
10.0 200.0 station_B 42
-10.0 -200.0 station_C
EOF

grt xy2geo -Ggrid_local.nc -Ogrid_geo.nc -C35.0/179.5
grt geo2xy -Ggrid_geo.nc -Ogrid_roundtrip.nc -C35.0/179.5
grt xy2geo -Gpoints_local.nc -Opoints_geo.nc -C-20.0/-179.5
grt geo2xy -Gpoints_geo.nc -Opoints_roundtrip.nc -C-20.0/-179.5
grt xy2geo -Gfinite_local.nc -Ofinite_geo.nc -C40.0/10.0
grt geo2xy -Gfinite_geo.nc -Ofinite_roundtrip.nc -C40.0/10.0
grt xy2geo -Qlocal_points.txt -Ogeo_points.txt -C35.0/179.5
grt geo2xy -Qgeo_points.txt -Olocal_points_roundtrip.txt -C35.0/179.5

if grt xy2geo -Ggrid_local.nc -Orejected_xy2geo.nc -R179.5/35.0 >/dev/null 2>&1; then
    echo "xy2geo unexpectedly accepted -R"
    exit 1
fi
if grt geo2xy -Ggrid_geo.nc -Orejected_geo2xy.nc -R179.5/35.0 >/dev/null 2>&1; then
    echo "geo2xy unexpectedly accepted -R"
    exit 1
fi

python -u test_xy2geo.py

rm -f grid_local.nc grid_geo.nc grid_roundtrip.nc \
    points_local.nc points_geo.nc points_roundtrip.nc \
    finite_local.nc finite_geo.nc finite_roundtrip.nc \
    local_points.txt geo_points.txt local_points_roundtrip.txt \
    rejected_xy2geo.nc rejected_geo2xy.nc

echo "test_xy2geo.sh: all checks passed"
