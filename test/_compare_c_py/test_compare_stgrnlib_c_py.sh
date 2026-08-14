#!/bin/bash
# Compare multi/single-depth STGRNLIB: CLI vs Python, and multi vs single slices

set -euo pipefail

modname="milrow"

# 粗网格以加快多深度测试
x1=-2; x2=2; dx=1
y1=-2; y2=2; dy=1

rm -rf stgrnlib_cmp
mkdir -p stgrnlib_cmp
cd stgrnlib_cmp

# (1) 单震源 + 单台站
grt static greenfn -M../../${modname} -D2/0 \
    -X$x1/$x2/$dx -Y$y1/$y2/$dy -e -Ostgrn_ss.nc

# (2) 多震源 + 单台站
grt static greenfn -M../../${modname} -Ds1,2,3 -Dr0 \
    -X$x1/$x2/$dx -Y$y1/$y2/$dy -e -Ostgrn_ms.nc

# (3) 单震源 + 多台站
grt static greenfn -M../../${modname} -Ds2 -Dr0,0.5,1 \
    -X$x1/$x2/$dx -Y$y1/$y2/$dy -e -Ostgrn_mr.nc

# (4) 多震源 + 多台站
grt static greenfn -M../../${modname} -Ds1,2,3 -Dr0,0.5,1 \
    -X$x1/$x2/$dx -Y$y1/$y2/$dy -e -Ostgrn_mm.nc

# 深度标签：整数写 0/1，小数点换成 p（与 compare_stgrnlib.py 一致）
depth_tag() {
    python -c "z=float('$1'); print(str(int(z)) if z==int(z) else str(z).replace('.','p'))"
}

# 用于多 vs 单逐层对比的参考单深度结果
for zs in 1 2 3; do
    for zr in 0 0.5 1; do
        zs_tag=$(depth_tag "$zs")
        zr_tag=$(depth_tag "$zr")
        grt static greenfn -M../../${modname} -D${zs}/${zr} \
            -X$x1/$x2/$dx -Y$y1/$y2/$dy -e \
            -Ostgrn_ref_zs${zs_tag}_zr${zr_tag}.nc
    done
done

cd -

python -u compare_stgrnlib.py

rm -rf stgrnlib_cmp
