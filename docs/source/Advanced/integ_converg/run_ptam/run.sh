#!/bin/bash
set -euo pipefail

rm -rf GRN* syn* stats* ptam* pygrtstats* static* *.png

# -------------------------------------------------------------------
# BEGIN DEPSRC 0.0 DGRN
# 使用 -Cp 指定使用 PTAM 进行收敛
grt greenfn -Mmilrow -D0/0 -N500/0.02 -OGRN -R5,8,10 -Cp -S50,100
grt ker2asc GRN_grtstats/milrow_0_0/K_0050_5.00000e+00 > stats
# 绘制图像部分见Python
# END DEPSRC 0.0 DGRN
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# BEGIN grt.ker2asc ptam
grt ker2asc GRN_grtstats/milrow_0_0/PTAM_0002_1.00000e+01/PTAM_0050_5.00000e+00 > ptam_stats
# END grt.ker2asc ptam
# -------------------------------------------------------------------

head -n 10 ptam_stats > ptam_stats_head
echo "..." >> ptam_stats_head







# -------------------------------------------------------------------
# BEGIN SGRN
# -S 表示输出核函数文件
grt static greenfn -Mmilrow -D0.1/0 -X2/2/1 -Y2/2/1 -Cp -S -Ostgrn.nc

# grt.ker2asc 也可以读取静态解输出的核函数文件，格式一致
grt ker2asc stgrtstats/milrow_0.1_0/K > static_stats

# 绘制图像部分见Python
# END SGRN
# -------------------------------------------------------------------

python run.py