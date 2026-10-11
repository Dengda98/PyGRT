#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt fault2nc -h

# 乱序的非方形网格，含零滑动和不同的行末延迟
grt fault2nc -Ffault_plane.inr -Ofault.nc
# 单块断层，支持 rake 和分量两种输入格式
grt fault2nc -Ffaults.inr -Ofault.nc
grt fault2nc -Ffault_single.inp -Ofault.nc

# 剖分模块生成的规则网格，可直接转换
grt faultdiv -Ffaults.inr -I0.4/0.6 -R2.5+o0.3/0.4 -Dc/1/1+d0.2 -Osubfaults.inp
grt fault2nc -Fsubfaults.inp -Ofault.nc

# 所有记录必须为 Kode=100
expect_fail grt fault2nc -Ffaults_rectangles.inp -Ofault.nc
# 网格不可缺块
sed '$d' fault_plane.inr > incomplete.inr
expect_fail grt fault2nc -Fincomplete.inr -Ofault.nc
# 输出不可覆盖输入
expect_fail grt fault2nc -Ffault_plane.inr -Ofault_plane.inr

rm -f fault.nc subfaults.inp incomplete.inr
remove_test_files
