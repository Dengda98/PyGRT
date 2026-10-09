#!/bin/bash

set -euo pipefail

expect_fail() {
    local desc="$1"
    shift
    set +e
    "$@" >/dev/null 2>&1
    local ret=$?
    set -e
    if [ "$ret" -eq 0 ]; then
        echo "ERROR: expected failure but succeeded: $desc" >&2
        exit 1
    fi
    echo "OK (failed as expected): $desc"
}

grt syn -h

grt greenfn -M../milrow -D2/3 -N40/0.02 -R10 -OGRN_SINGLE -s
grt syn -GGRN_SINGLE -A22 -S1e20 -Osyn_single
grt syn -GGRN_SINGLE -Ds2 -Dr3 -R10 -A22 -S1e20 -Osyn_single_explicit
expect_fail "wrong single source depth is rejected" \
    grt syn -GGRN_SINGLE -Ds2.1 -A22 -S1e20 -Osyn_bad
expect_fail "wrong single receiver depth is rejected" \
    grt syn -GGRN_SINGLE -Dr3.1 -A22 -S1e20 -Osyn_bad
expect_fail "wrong single epicentral distance is rejected" \
    grt syn -GGRN_SINGLE -R11 -A22 -S1e20 -Osyn_bad

grt greenfn -M../milrow -Ds1,2 -Dr3 -N40/0.02 -R10 -OGRN_SRC_MULTI -s
grt syn -GGRN_SRC_MULTI -Ds1 -A22 -S1e20 -Osyn_src_multi
grt syn -GGRN_SRC_MULTI -Ds1 -Dr3 -A22 -S1e20 -Osyn_src_multi_explicit
expect_fail "multiple source depths require -Ds" \
    grt syn -GGRN_SRC_MULTI -A22 -S1e20 -Osyn_bad
expect_fail "wrong single receiver depth is rejected in source-depth library" \
    grt syn -GGRN_SRC_MULTI -Ds1 -Dr3.1 -A22 -S1e20 -Osyn_bad

grt greenfn -M../milrow -Ds2 -Dr1,3 -N40/0.02 -R10 -OGRN_RCV_MULTI -s
grt syn -GGRN_RCV_MULTI -Dr1 -A22 -S1e20 -Osyn_rcv_multi
grt syn -GGRN_RCV_MULTI -Ds2 -Dr1 -A22 -S1e20 -Osyn_rcv_multi_explicit
expect_fail "multiple receiver depths require -Dr" \
    grt syn -GGRN_RCV_MULTI -A22 -S1e20 -Osyn_bad
expect_fail "wrong single source depth is rejected in receiver-depth library" \
    grt syn -GGRN_RCV_MULTI -Dr1 -Ds2.1 -A22 -S1e20 -Osyn_bad

grt greenfn -M../milrow -D2/3 -N40/0.02 -R5,10 -OGRN_DIST_MULTI -s
grt syn -GGRN_DIST_MULTI -R10 -A22 -S1e20 -Osyn_dist_multi
expect_fail "multiple epicentral distances require -R" \
    grt syn -GGRN_DIST_MULTI -A22 -S1e20 -Osyn_bad

# 分次计算时补齐三个采样轴的所有组合，测试 syn 的根目录查询
grt greenfn -M../milrow -Ds1,2 -Dr0,3 -N600/0.02 -R5,10 -OGRN -s
grt greenfn -M../milrow -D2/3 -N600/0.02 -R10 -e -OGRN

grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Osyn 
grt syn -GGRN -Ds1 -Dr0 -R10 -A22 -S1e20 -Osyn_multi_1_0_10
grt syn -GGRN -Ds2 -Dr3 -R10 -A22 -S1e20 -Osyn_multi_2_3_10
expect_fail "selectors cannot be used with a subdirectory" \
    grt syn -GGRN/milrow_2_3_10 -Ds2 -Dr3 -R10 -A22 -S1e20 -Osyn_bad
grt syn -GGRN -Ds2 -Dr3 -R9 -A22 -S1e20 -Osyn_root_linear
grt syn -GGRN -Ds2 -Dr3 -R9 -A22 -S1e20 -i0 -Osyn_root_nearest
expect_fail "root queries reject distances outside the library" \
    grt syn -GGRN -Ds2 -Dr3 -R11 -A22 -S1e20 -i0 -Osyn_bad
for option in -Ds2 -Dr3 -R10; do
    expect_fail "node directory rejects $option independently" \
        grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 "$option" -Osyn_bad
done
for mode in 0 1; do
    grt syn -GGRN/milrow_2_3_10/ -A22 -S1e20 -i"$mode" -Osyn_node_i"$mode"
done
expect_fail "root directory requires -R when depths match multiple distances" \
    grt syn -GGRN -Ds2 -Dr3 -A22 -S1e20 -Osyn_bad
grt syn -GGRN/milrow_2_3_10 -A22 -S1e16 -F-1/2/-4 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -M33/44/55 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -Su1e10 -M33/44/55 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -Su1e10 -M33/44 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20  -T1/-2/-5/0.5/3/1.2 -Osyn 

# 震源时间函数按矩形法进行面积归一化
# 自定义时间函数检查 dt 乘样本和，非单位面积时警告并自动归一化
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Dp/0.6 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Dt/0.2/0.2/0.3 -Osyn
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Dt/0.4/0/0.4 -Osyn
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Dc/0.2/0.4 -Osyn
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -Dr/1.2 -Osyn
cat > tfile <<EOF
0.0
5.0
10.0
5.0
10.0
10.0
10.0
0.0
EOF
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -D0/tfile -Osyn_custom
cat > tfile_invalid <<EOF
0.0 0.1
EOF
expect_fail "custom time function file must contain exactly one column" \
    grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -D0/tfile_invalid -Osyn_custom_invalid
cat > tfile_warning <<EOF
-0.2
0.6
EOF
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -D0/tfile_warning -Osyn_custom_warning > custom_warning.log
grep -q "Custom time function sequence sum" custom_warning.log
rm -f tfile tfile_invalid tfile_warning custom_warning.log

grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -M33/44/55 -I1 -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -M33/44/55 -J1 -Osyn 

grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -N -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -e -Osyn 
grt syn -GGRN/milrow_2_3_10 -A22 -S1e20 -N -e -Osyn 

# Kode=100 的 rake 接收断层在零滑动量时仍保留滑动角，省略剖分时取中心点
cat > rcv_faults.inr <<'EOF'
# X-start Y-start X-fin Y-fin Kode rake netslip dip top bot
xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx
1 8 0 10 0 100 45 0 90 1 3
2 8 0 10 0 100 -30 1 90 1 3
3 8 0 10 0 100 90 0 90 1 3
EOF
grt syn -GGRN -Ds2 -S1e20 -Urcv_faults.inr -Osyn_rcv_center
grt syn -GGRN -Ds2 -S1e20 -Urcv_faults.inr+i1/1 -Osyn_rcv_subdiv
# 单点接收文件也使用子目录
head -n 3 rcv_faults.inr > rcv_fault_single.inr
printf '0 10 3\n' > rcv_single.txt
grt syn -GGRN -Ds2 -S1e20 -Urcv_fault_single.inr -Osyn_rcv_fault_single
grt syn -GGRN -Ds2 -S1e20 -Qrcv_single.txt -Osyn_rcv_point_single
expect_fail "receiver subdivision sizes must be positive" \
    grt syn -GGRN -Ds2 -S1e20 -Urcv_faults.inr+i0/1 -Osyn_bad
rm -f rcv_faults.inr rcv_fault_single.inr rcv_single.txt
rm -rf syn_rcv_center syn_rcv_subdiv syn_rcv_fault_single syn_rcv_point_single


python -u test_syn.py


rm -rf GRN GRN_SINGLE GRN_SRC_MULTI GRN_RCV_MULTI GRN_DIST_MULTI \
    syn syn_custom syn_custom_warning syn_custom_invalid syn_single syn_bad syn_src_multi syn_rcv_multi syn_dist_multi \
    syn_single_explicit syn_src_multi_explicit syn_rcv_multi_explicit \
    syn_subdir_bad syn_multi_1_0_10 syn_multi_2_3_10 syn_root_linear syn_root_nearest syn_node_i0 syn_node_i1
