#!/bin/bash
set -euo pipefail

# 构建有限震源输入
cat > fault.inp << 'EOF'
  #  X-start  Y-start    X-fin    Y-fin Kode   rt.lat  reverse      dip      top      bot
--- -------- -------- -------- -------- ---- -------- -------- -------- -------- --------
  1   -4.000    0.000    4.000    0.000  100    2.000    0.000   60.000    2.000    5.464
EOF

# BEGIN GRN
# 接收点位于地表，-e 增加位移空间偏导数
grt static greenfn -Mmod3 -Ds2/6/1 -Dr0 -R0/36/2 -Ostgrn.nc -e
# END GRN

# BEGIN SRC SYN
# 使用有限震源，地表二维网格输出 ZNE 位移及偏导数
grt static syn -Gstgrn.nc -Cfault.inp+i1/1 -Dr0 -N -X-20/20/1 -Y-20/20/1 -Ostsyn_ff.nc -e
# END SRC SYN

# BEGIN COULOMB
grt static stress -Gstsyn_ff.nc
grt static sproj -Gstsyn_ff.nc -M90/60/0
grt static coulomb -Gstsyn_ff.nc -F0.4
# END COULOMB

# 绘制结果图
python plot.py
