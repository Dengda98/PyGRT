#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt greenfn -h

# 采样、频带、积分控制和时间偏移
grt greenfn -M../milrow -D2/0 -N32/0.1 -R5 -OGRN
grt greenfn -M../milrow -D2/0 -N32/0.1+w0.6+n2+a+f -R5 -e -P2 -OGRN
grt greenfn -M../milrow -D2/0 -N32/0.1 -R5 -H0.5/2 -L10 -OGRN
grt greenfn -M../milrow -D2/0 -N32/0.1 -R5 -K+k4+f+s1.2+e1e-3+v1.5 -OGRN
grt greenfn -M../milrow -D2/0 -N32/1 -R100 -L+l10+o2 -E-2/9 -OGRN
grt greenfn -M../milrow -D2/0 -N32/1 -R100 -L+a1e-3 -Ep-1 -OGRN

# 收敛方法、边界条件、源类型和积分统计
grt greenfn -M../milrow -D2/0 -N32/0.1 -R5 -Cd -BrF -OGRN
grt greenfn -M../milrow -D2/0 -N32/0.1 -R5 -Cp -BhR -OGRN
grt greenfn -M../milrow -D2/0 -N32/0.1 -R5 -Cn -BrH -Gev -S1,2 -OGRN

# 距离和深度的列表、范围、文件输入
grt greenfn -M../milrow -D2/0 -N16/0.1 -R2,5,10 -OGRN
grt greenfn -M../milrow -Ds1/2/1 -Dr0,1 -N16/0.1 -R2/10/4 -OGRN
grt greenfn -M../milrow -Dsdepsrc -Drdeprcv -N16/0.1 -Rdists -OGRN -s

# 震中距列表不可包含逆序或重复值
expect_fail grt greenfn -M../milrow -D2/0 -N32/0.1 -R5,2 -OGRN
# 不可同时设置单深度和多深度参数
expect_fail grt greenfn -M../milrow -D2/0 -Ds1,2 -Dr0 -N32/0.1 -R5 -OGRN
# 未指定多深度参数对应的接收深度
expect_fail grt greenfn -M../milrow -Ds1,2 -N32/0.1 -R5 -OGRN
# 时间采样间隔不可为零或负数
expect_fail grt greenfn -M../milrow -D2/0 -N32/0 -R5 -OGRN

# 频带下限不可超过上限
expect_fail grt greenfn -M../milrow -D2/0 -N32/0.1 -R5 -H2/0.5 -OGRN
# 不可同时设置 FIM 和 SAFIM
expect_fail grt greenfn -M../milrow -D2/0 -N32/0.1 -R5 -L+l10+a1e-3 -OGRN

python -u test_greenfn.py

rm -rf GRN GRN_grtstats
remove_test_files
