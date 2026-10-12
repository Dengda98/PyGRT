#!/bin/bash
set -euo pipefail

python prepare_points.py
head -n 10 rcv_pts.txt > rcv_head
echo "..." >> rcv_head

# BEGIN GRN
grt static greenfn -Mmilrow -D2/0 -R0/10/0.1 -Ostgrn.nc
# END GRN

# BEGIN SYN GRID
# 爆炸源
grt static syn -Gstgrn.nc -S1e24 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_ex.nc -N

# 单力源
grt static syn -Gstgrn.nc -S1e16 -F1/-0.5/2 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_sf.nc -N

# 剪切源
grt static syn -Gstgrn.nc -S1e24 -M33/50/120 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_dc.nc -N

# 剪切源
grt static syn -Gstgrn.nc -S1e24 -M33/90/0 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_dc2.nc -N

# 张裂源
grt static syn -Gstgrn.nc -S1e24 -M33/50 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_ts.nc -N

# 张裂源
grt static syn -Gstgrn.nc -S1e24 -M33/90 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_ts2.nc -N

# 矩张量源
grt static syn -Gstgrn.nc -S1e24 -T0.1/-0.2/1.0/0.3/-0.5/-2.0 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_mt.nc -N

# 矩张量源
grt static syn -Gstgrn.nc -S1e24 -T0/-0.2/0/0/0/0 -X-3/3/0.15 -Y-2.5/2.5/0.15 -Ostsyn_mt2.nc -N
# END SYN GRID

# BEGIN SYN POINTS
# 使用 -Q 来传入任意点坐标文件
grt static syn -Gstgrn.nc -S1e24 -M33/90/0 -Qrcv_pts.txt -Ostsyn_points.nc -N
# END SYN POINTS

python run.py
bash plot.sh
python plot.py

# 打包本页文档、输入文件、脚本和运行结果
rm -f *.tar.gz
example=$(basename "$PWD")
tar -czf "../${example}.tar.gz" -C .. "$example"
mv "../${example}.tar.gz" "${example}.tar.gz"

# 打包后删除计算结果，节省本地空间
rm -f stgrn.nc stsyn_*.nc rcv_pts.txt syn_grid.pdf
