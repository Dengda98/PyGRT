#!/bin/bash 

set -euo pipefail

rm -rf GRN* pygrtstats* *.png

# -----------------------------------------------------------------
# BEGIN GRN
# -S 后不指定索引表示输出所有频率点的核函数
# -Cn 禁用收敛算法
# -L20 定义波数积分间隔dk
grt greenfn -Mmod1 -D0.03/0 -N500/0.02 -OGRN -R1 -K+v0.1 -S -L20 -Cn
# END GRN
# -----------------------------------------------------------------

python run.py

rm -rf GRN* pygrtstats*