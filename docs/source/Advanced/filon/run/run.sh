#!/bin/bash
# BEGIN
rm -rf GRN*

# method 1: DWM
# 使用-H来节省一些时间，只输出单个频率以作示例
grt greenfn -Mmilrow -OGRN1 -N2000/1 -D5/0 -R2500 -L20 -H0.25/0.25  -S

# method 2: SAFIM
# 使用-La指定相关参数
grt greenfn -Mmilrow -OGRN2 -N2000/1 -D5/0 -R2500 -La20/1e-2/10 -S500
# END

python run.py
python plot.py