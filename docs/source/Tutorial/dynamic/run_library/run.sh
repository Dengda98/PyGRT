#!/bin/bash
set -euo pipefail

# 构建有限震源输入
cat > fault.inp << 'EOF'
  #  X-start  Y-start    X-fin    Y-fin Kode   rt.lat  reverse      dip      top      bot
--- -------- -------- -------- -------- ---- -------- -------- -------- -------- --------
  1   -4.000    0.000    4.000    0.000  100    2.000    0.000   60.000    2.000    5.464
EOF

# 生成子断层、时间函数和从断层面中心向外传播的破裂延迟
grt faultdiv -Ffault.inp -I1/1 -R2.5+o4/2 -Dc/1/1 -Ofault_rupture.inp

# 重新建库，避免保留旧的接收深度采样节点
rm -rf GRN

# BEGIN GRN
# -N 的 +a 表示计算所有频点
grt greenfn -Mmod3 -Ds2/6/1 -Dr0 -R2/36/2 -N512/0.1+a -L+a1e-3 -OGRN -e
# END GRN

# BEGIN SRC SYN
# 每行已经是 1×1 km 子断层，+i0/0 保证每行只取一个子源中心
# 不设置全局 -D，以保留文件中的破裂延迟
grt syn -GGRN -Cfault_rupture.inp+i0/0 -R20 -A30 -Dr0 -Osyn_ff -e -I1
grt stress -Gsyn_ff
grt sproj -Gsyn_ff -M90/60/0
grt coulomb -Gsyn_ff -F0.4
# END SRC SYN

python plot.py
