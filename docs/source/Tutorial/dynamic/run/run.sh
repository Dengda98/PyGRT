#!/bin/bash 

# BEGIN GRN
# 不同震源深度、接收点深度和震中距的格林函数会在 
# GRN/milrow_{depsrc}_{deprcv}_{dist}/ 路径下，使用SAC格式保存。 
grt -Mmilrow -D2/0 -N500/0.02 -OGRN -R5,8,10
# END GRN

# BEGIN grt.b2a
grt.b2a GRN/milrow_2_0_10/HFZ.sac > HFZ
# END grt.b2a
head -n 10 HFZ > HFZ_head
echo "..." >> HFZ_head

# BEGIN SYN EXP
# 合成结果在 syn_exp/ 目录下，以SAC格式保存。
grt.syn -GGRN/milrow_2_0_10 -S1e24 -A30 -Osyn_exp
# END SYN EXP

# BEGIN SYN SF
# 合成结果在 syn_sf/ 目录下，以SAC格式保存。
grt.syn -GGRN/milrow_2_0_10 -S1e16 -A30 -F1/-0.5/2 -Osyn_sf
# END SYN SF

# BEGIN SYN DC
# 合成结果在 syn_dc/ 目录下，以SAC格式保存。
grt.syn -GGRN/milrow_2_0_10 -S1e24 -A30 -M33/50/120 -Osyn_dc
# END SYN DC

# BEGIN SYN MT
# 合成结果在 syn_mt/ 目录下，以SAC格式保存。
grt.syn -GGRN/milrow_2_0_10 -S1e24 -A30 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -Osyn_mt
# END SYN MT


# BEGIN ZNE
# 传入-N可输出ZNE分量
grt.syn -GGRN/milrow_2_0_10 -S1e24 -A30 -M33/50/120 -Osyn_dc_zne -N
# END ZNE


# BEGIN TIME FUNC
# -Dt定义梯形波的三个时间参数，爬升截止时刻/平稳截止时刻/结束时刻
# 前两个时刻重合则退化为三角波
grt.syn -GGRN/milrow_2_0_10 -S1e16 -A30 -F1/-0.5/2 -Osyn_sf_trig -Dt/0.3/0.3/0.6
# END TIME FUNC


# BEGIN INT DIF
# 积分，1表示积分1次
grt.syn -GGRN/milrow_2_0_10 -S1e24 -A30 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -I1 -Osyn_mt_int1
# 微分，1表示微分1次
grt.syn -GGRN/milrow_2_0_10 -S1e24 -A30 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -J1 -Osyn_mt_dif1
# END INT DIF