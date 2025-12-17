#!/bin/bash
set -euo pipefail

rm -rf GRN* syn* stats* pygrtstats* *.png

# -------------------------------------------------------------------
# BEGIN DGRN
# -S<i1>,<i2>,...表示输出对应频率点索引值的核函数文件
grt greenfn -Mmilrow -D2/0 -N500/0.02 -OGRN -R5,8,10 -S50,100
# END DGRN
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# BEGIN grt.ker2asc
grt ker2asc GRN_grtstats/milrow_2_0/K_0050_5.00000e+00 > stats
# END grt.ker2asc
# -------------------------------------------------------------------

head -n 10 stats > stats_head
echo "..." >> stats_head


# -------------------------------------------------------------------
# BEGIN DEPSRC 0.0 DGRN
# 使用 -Cn 禁用任何收敛算法
grt greenfn -Mmilrow -D0/0 -N500/0.02 -OGRN -R5,8,10 -Cn -S50,100
grt ker2asc GRN_grtstats/milrow_0_0/K_0050_5.00000e+00 > stats
# 绘制图像部分见Python
# END DEPSRC 0.0 DGRN
# -------------------------------------------------------------------

python run.py

rm -rf GRN* syn* pygrtstats*