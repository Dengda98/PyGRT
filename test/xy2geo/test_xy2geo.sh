#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt xy2geo -h
grt geo2xy -h

# 文本坐标
grt xy2geo -Qrcv_points.txt -Ogeo.txt -C35/179.5
grt geo2xy -Qgeo.txt -Olocal.txt -C35/179.5

# 用普通模块准备网格、任意接收点和接收断层文件
grt okada -I6/3.464/2.7 -S1e20 -Ds10 -Dr0 -X-2/2/2 -Y-2/2/2 -Ogrid.nc
grt xy2geo -Ggrid.nc -Ogeo.nc -C35/179.5
grt geo2xy -Ggeo.nc -Olocal.nc -C35/179.5
grt okada -I6/3.464/2.7 -S1e20 -Ds10 -Qrcv_points.txt -Opoints.nc
grt xy2geo -Gpoints.nc -Ogeo.nc -C-20/-179.5
grt geo2xy -Ggeo.nc -Olocal.nc -C-20/-179.5
grt okada -I6/3.464/2.7 -S1e20 -Ds10 -Urcv_faults.inr -Ofaults.nc
grt xy2geo -Gfaults.nc -Ogeo.nc -C40/10
grt geo2xy -Ggeo.nc -Olocal.nc -C40/10

# 不可同时设置 NetCDF 文件和文本坐标文件
expect_fail grt xy2geo -Ggrid.nc -Qrcv_points.txt -Ogeo.nc -C35/10
# 参考纬度不可达到或超过 ±90 度
expect_fail grt xy2geo -Qrcv_points.txt -Ogeo.txt -C91/10
# 未指定参考经纬度
expect_fail grt geo2xy -Qgeo.txt -Olocal.txt

# 参考经度不可超出 [-180, 180] 度
expect_fail grt xy2geo -Qrcv_points.txt -Ogeo.txt -C35/181

python -u test_xy2geo.py

rm -rf grid.nc points.nc faults.nc geo.nc local.nc geo.txt local.txt
remove_test_files
