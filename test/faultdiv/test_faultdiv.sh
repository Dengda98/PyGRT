#!/bin/bash

set -euo pipefail
source ../common.sh
create_test_files

grt faultdiv -h

# 分量格式和 rake 格式的矩形断层，尺寸调整为整数等分
grt faultdiv -Ffaults_rectangles.inp -I0.4/0.6 -Osubfaults_rectangles.inp
grt faultdiv -Ffaults.inr -I0.7/0.8 -Osubfaults.inp
# 输入行末的时间函数和延迟忽略
grt faultdiv -Ffaults_negative_delay.inr -I0.4/0.6 -Osubfaults.inp

# 径向破裂和平面破裂，可叠加整体延迟
grt faultdiv -Ffaults.inr -I0.4/0.6 -R2.5+o0.3/0.4 -Dc/0.2/0.3+d0.2 -Osubfaults.inp
grt faultdiv -Ffaults_rectangles.inp -I0.4/0.6 -R2.5+d90 -Osubfaults_rectangles.inp
# 仅设置时间函数，所有子断层采用相同延迟
grt faultdiv -Ffaults.inr -I0.4/0.6 -Dp/0.5+d1.2 -Osubfaults.inp

# 地表倾斜断层的剖分结果可直接合成，保留已有剖分和行末时间函数
grt faultdiv -Ffaults_surface.inr -I5/5 -R3.5+d90 -Dt/0/0.5/0 -Osubfaults_surface.inp
grt lamb -H8/4.62/3.3 -N32/0.1 -Csubfaults_surface.inp+i0/0 -R20 -A0 -Dr2 -I1 -P2 -Osyn_surface

# 不支持点源
expect_fail grt faultdiv -Ffaults.inp -I1/1 -Osubfaults.inp
# 剖分尺寸须为正数
expect_fail grt faultdiv -Ffaults.inr -I0/1 -Osubfaults.inp
# 径向破裂要求输入只有一条断层记录
expect_fail grt faultdiv -Ffaults_rectangles.inp -I1/1 -R2.5+o0/0 -Osubfaults.inp
# 起裂点须在断层面内
expect_fail grt faultdiv -Ffaults.inr -I1/1 -R2.5+o2/0 -Osubfaults.inp
# 输出不可覆盖输入
expect_fail grt faultdiv -Ffaults.inr -I1/1 -Ofaults.inr

rm -f subfaults_rectangles.inp subfaults.inp subfaults_surface.inp
rm -rf syn_surface
remove_test_files
