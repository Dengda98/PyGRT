#!/bin/bash

set -euo pipefail

# 2×2 库仅包含四个角点，目标位置的精确格林函数单独计算
# 两次计算使用相同的 dk = 2π/240 km，避免积分采样差异影响比较
grt greenfn -Mmilrow -Ds4.5,5.5 -Dr0 -R9,11 -N1000/0.02 -OGRN -s
grt greenfn -Mmilrow -D5/0 -R10 -N1000/0.02 -OGRN_EXACT -s

# 剪切源卷积梯形时间函数，再积分一次得到阶跃响应，输出 ZRT 分量
SET="-A30 -S1e23 -M33/50/120 -Dt/0.2/0.2/0.2 -I1"
grt syn -GGRN -Ds5 -Dr0 -R10 $SET -Osyn_p_aligned
grt syn -GGRN -Ds5 -Dr0 -R10 $SET -i0 -Osyn_nearest
grt syn -GGRN_EXACT/milrow_5_0_10 $SET -Osyn_exact

# 普通线性插值的四个角点，各自使用相同的源机制和时间处理
mkdir -p syn_linear
for depth in 4.5 5.5; do
    for distance in 9 11; do
        grt syn -GGRN/milrow_${depth}_0_${distance} $SET -Osyn_linear/${depth}_${distance}
    done
done

python plot.py

cp comparison.svg cover.svg

ex=$(basename "$(pwd)")
rm -f "$ex.tar.gz"
cd .. && tar -czvf "${ex}.tar.gz" "${ex}" && mv "${ex}.tar.gz" "${ex}" && cd -

rm -rf GRN GRN_EXACT syn_p_aligned syn_nearest syn_exact syn_linear
