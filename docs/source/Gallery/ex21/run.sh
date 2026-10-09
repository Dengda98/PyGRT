#!/bin/bash

set -euo pipefail

depsrc=3
deprcv=0
dist=10
azimuth=60
strike=33
dip=50
rake=120
moment=1e24
nt=2500
dt=0.01

# 静态解用 -X/-Y 指定单点，北向和东向坐标由距离、方位角换算，单位为 km
read -r north east < <(awk -v r="$dist" -v a="$azimuth" \
    'BEGIN { a *= atan2(0, -1) / 180; printf "%.15g %.15g\n", r * cos(a), r * sin(a) }')

for model in layered halfspace; do
    # 保留全部低频，避免时间积分后的永久位移发生漂移
    grt greenfn -M$model -D$depsrc/$deprcv -R$dist -N$nt/$dt+a -Gs -OGRN_$model
    grt syn -GGRN_$model -Ds$depsrc -Dr$deprcv -R$dist -A$azimuth -S$moment -M$strike/$dip/$rake \
        -Dt/0.1/0.1/0.1 -I1 -N -Osyn_$model

    grt static greenfn -M$model -D$depsrc/$deprcv -R$dist -L60 -Ostatic_greenfn_$model.nc
    grt static syn -Gstatic_greenfn_$model.nc -S$moment -M$strike/$dip/$rake -Dr$deprcv \
        -X$north/$north/1 -Y$east/$east/1 -N -Ostatic_$model.nc
done

# 均匀半空间的动态、静态解析解，震源和台站参数与数值解一致
grt lamb -H6/3.464/2.7 -Ds$depsrc -Dr$deprcv -R$dist -A$azimuth -N$nt/$dt \
    -S$moment -M$strike/$dip/$rake -Dt/0.1/0.1/0.1 -I1 -n -Olamb
grt okada -I6/3.464/2.7 -Ds$depsrc -Dr$deprcv -X$north/$north/1 -Y$east/$east/1 \
    -S$moment -M$strike/$dip/$rake -N -Ookada.nc

python plot.py

cp compare_layered.svg cover.svg

rm -rf GRN_layered GRN_halfspace syn_layered syn_halfspace lamb
rm -f static_greenfn_layered.nc static_greenfn_halfspace.nc static_layered.nc static_halfspace.nc okada.nc

ex=$(basename "$(pwd)")
rm -f "$ex.tar.gz"
cd .. && tar -czvf "$ex.tar.gz" "$ex" && mv "$ex.tar.gz" "$ex" && cd -
