#!/bin/bash
set -euo pipefail

# BEGIN GRN
# 不同震源深度、接收点深度和震中距的格林函数会在
# GRN/milrow_{depsrc}_{deprcv}_{dist}/ 路径下，使用SAC格式保存
grt greenfn -Mmilrow -D2/0 -N500/0.02 -OGRN -R5,8,10
# END GRN

# BEGIN SYN SOURCES
# 合成结果在对应目录下，以SAC格式保存

# 爆炸源，标量矩 1e24 dyne·cm
grt syn -GGRN/milrow_2_0_10 -S1e24 -A30 -Osyn_ex
# 也可以将 -G 指向库根目录，按目标几何查询，默认线性插值，-i0 改为最近邻
# 本库的源深和台深各只有一个值，-Ds/-Dr 可省略，-R7 在 5 km 和 8 km 节点之间插值
# grt syn -GGRN -R7 -S1e24 -A30 -Osyn_ex

# 单力源，(fN, fE, fZ)=(1, -0.5, 2)【北向，东向，垂直向下】，标度 1e16 dyne
grt syn -GGRN/milrow_2_0_10 -S1e16 -A30 -F1/-0.5/2 -Osyn_sf

# 剪切源，strike=33°，dip=50°，rake=120°，标量矩 1e24 dyne·cm
grt syn -GGRN/milrow_2_0_10 -S1e24 -A30 -M33/50/120 -Osyn_dc

# 走滑剪切源，strike=33°，dip=90°，rake=0°，标量矩 1e24 dyne·cm
grt syn -GGRN/milrow_2_0_10 -S1e24 -A30 -M33/90/0 -Osyn_dc2

# 张裂源，strike=33°，dip=50°，标量矩 1e24 dyne·cm
grt syn -GGRN/milrow_2_0_10 -S1e24 -A30 -M33/50 -Osyn_ts

# 矩张量源，标量矩 1e24 dyne·cm， -T 参数的每一项的下标中，X为北向，Y为东向，Z为垂直向下
grt syn -GGRN/milrow_2_0_10 -S1e24 -A30 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -Osyn_mt
# END SYN SOURCES

# BEGIN ZNE
# 传入-N可输出ZNE分量
grt syn -GGRN/milrow_2_0_10 -S1e24 -A30 -M33/50/120 -Osyn_dc_zne -N
# END ZNE

# BEGIN TIME FUNC
# -Dt定义梯形波的三个时间参数，上升段时长/平台段时长/下降段时长
# 平台段时长为零则退化为三角波
grt syn -GGRN/milrow_2_0_10 -S1e16 -A30 -F1/-0.5/2 -Osyn_sf_trig -Dt/0.3/0/0.3
# -Dc定义非对称余弦波，上升段时长/下降段时长
grt syn -GGRN/milrow_2_0_10 -S1e16 -A30 -F1/-0.5/2 -Osyn_sf_cosine -Dc/0.1/0.5
# END TIME FUNC

# BEGIN INT DIF
# 积分，1表示积分1次
grt syn -GGRN/milrow_2_0_10 -S1e24 -A30 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -I1 -Osyn_mt_int1
# 微分，1表示微分1次
grt syn -GGRN/milrow_2_0_10 -S1e24 -A30 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -J1 -Osyn_mt_dif1
# END INT DIF

python run.py
python plot.py

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -rf GRN syn_ex syn_sf syn_dc syn_dc2 syn_ts syn_mt syn_dc_zne syn_sf_trig syn_sf_cosine syn_mt_int1 syn_mt_dif1
