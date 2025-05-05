#!/bin/bash
set -eu

rm -rf GRN* syn* stats* ptam* pygrtstats* static* *.png

# -------------------------------------------------------------------
# BEGIN DGRN
# -S<i1>,<i2>,...表示输出对应频率点索引值的核函数文件
grt -Mmilrow -D2/0 -N500/0.02 -OGRN -R5,8,10 -S50,100
# END DGRN
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# BEGIN grt.k2a
grt.k2a GRN_grtstats/milrow_2_0/K_0050_5.00000e+00 > stats
# END grt.k2a
# -------------------------------------------------------------------

head -n 10 stats > stats_head
echo "..." >> stats_head


# -------------------------------------------------------------------
# BEGIN DEPSRC 0.0 DGRN
grt -Mmilrow -D0/0 -N500/0.02 -OGRN -R5,8,10 -S50,100
grt.k2a GRN_grtstats/milrow_0_0/K_0050_5.00000e+00 > stats
# 绘制图像部分见Python
# END DEPSRC 0.0 DGRN
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# BEGIN grt.k2a ptam
grt.k2a GRN_grtstats/milrow_0_0/PTAM_0002_1.00000e+01/PTAM_0050_5.00000e+00 > ptam_stats
# END grt.k2a ptam
# -------------------------------------------------------------------

head -n 10 ptam_stats > ptam_stats_head
echo "..." >> ptam_stats_head







# -------------------------------------------------------------------
# BEGIN SGRN
# -S 表示输出核函数文件
stgrt -Mmilrow -D0.1/0 -X2/2/1 -Y2/2/1 -S > grn

# grt.k2a 也可以读取静态解输出的核函数文件，格式一致
grt.k2a stgrtstats/milrow_0.1_0/K > static_stats

# 绘制图像部分见Python
# END SGRN
# -------------------------------------------------------------------
