#!/bin/bash

set -euo pipefail

XY="-X-200/200/10 -Y-210/210/10"

# 参考点经纬度
lat0="-8.3101000"
lon0="121.3517000"

# 分层模型中的解
OUT="stsyn_ff.nc"
grt static syn -Gstgrn.nc -Ccoulomb-fault.inp -Dr0 $XY -e -N -O$OUT

grt static stress $OUT
grt static sproj -G$OUT -M90/32/90
grt static coulomb -G$OUT -F0.4
grt xy2geo -G$OUT -C$lat0/$lon0 -O${OUT}.geo

# Okada 解
OUT="okada_ff.nc"
grt okada -I6/3.464/2.7 -Ccoulomb-fault.inp -Dr0 $XY -e -N -O$OUT

grt static stress $OUT
grt static sproj -G$OUT -M90/32/90
grt static coulomb -G$OUT -F0.4
grt xy2geo -G$OUT -C$lat0/$lon0 -O${OUT}.geo

