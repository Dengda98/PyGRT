#!/bin/bash


# BEGIN GRN
# -e 表示计算空间导数
grt -Mmilrow -D2/0 -N500/0.02 -OGRN -R10 -e
# END GRN

# BEGIN SYN DC
# -e 表示计算空间导数
grt.syn -GGRN/milrow_2_0_10 -S1e24 -A30 -M33/50/120 -Osyn_dc -e
# END SYN DC

# BEGIN ZNE
# -N 表示输出ZNE分量
grt.syn -GGRN/milrow_2_0_10 -S1e24 -A30 -M33/50/120 -Osyn_dc_zne -e -N
# END ZNE


# BEGIN STRAIN 
# 指定文件夹以及中间名，会在原文件夹内输出SAC格式的应变
grt.strain syn_dc_zne/out
# END STRAIN

# BEGIN STRESS 
# 指定文件夹以及中间名，会在原文件夹内输出SAC格式的应变
grt.stress syn_dc_zne/out
# END STRESS