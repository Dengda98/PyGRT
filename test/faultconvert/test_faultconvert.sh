#!/bin/bash

set -euo pipefail

source ../common.sh

# 每次从 USGS 下载同一产品版本的两种格式，运行后删除
url="https://earthquake.usgs.gov/product/finite-fault/us6000tkt2_1/us/1787347333870"
for file in basic_inversion.param complete_inversion.fsp; do
    curl -fLsS --retry 3 "$url/$file" -o "$file"
done

# 历史多段模型使用较低精度的边界坐标，两种格式均保持按段导出的顺序
url="https://earthquake.usgs.gov/product/finite-fault/us1000778i/us/1539813034137"
curl -fLsS --retry 3 "$url/basic_inversion.param"  -o legacy.param
curl -fLsS --retry 3 "$url/complete_inversion.fsp" -o legacy.fsp

grt faultconvert -h

grt faultconvert -Fp/basic_inversion.param -C-8.3101/121.3517 -Oconverted_param.inr
grt faultconvert -Ff/complete_inversion.fsp -C-8.3101/121.3517 -Oconverted_fsp.inr
grt faultconvert -Fp/basic_inversion.param -C-8.3101/121.3517 -N1 -Oconverted_param.inr

# 两种格式使用同一首块中心，均可选择不同的局部坐标原点
grt faultconvert -Fp/basic_inversion.param -C-8/122 -Oreference_param.inr
grt faultconvert -Ff/complete_inversion.fsp -C-8/122 -Oreference_fsp.inr

# 实际记录数与网格一致时，Nsbfs 声明偏小或偏大均可警告后转换
sed -E 's/Nsbfs[[:space:]]*=[[:space:]]*[0-9]+/Nsbfs = 1/' complete_inversion.fsp > count.fsp
grt faultconvert -Ff/count.fsp -C-8.3101/121.3517 -Oconverted_fsp.inr
sed -E 's/Nsbfs[[:space:]]*=[[:space:]]*[0-9]+/Nsbfs = 999/' complete_inversion.fsp > count.fsp
grt faultconvert -Ff/count.fsp -C-8.3101/121.3517 -Oconverted_fsp.inr

# 历史坐标不需要满足投影后的边界残差限制，多段可整体或选段转换
grt faultconvert -Fp/legacy.param -C-42.757/173.077 -Olegacy.inr
grt faultconvert -Ff/legacy.fsp   -C-42.757/173.077 -Olegacy.inr
grt faultconvert -Fp/legacy.param -C-42.757/173.077 -N3 -Olegacy.inr
grt fault2nc -Flegacy.inr -Ofault.nc
grt faultconvert -Ff/legacy.fsp   -C-42.757/173.077 -N3 -Olegacy.inr
grt fault2nc -Flegacy.inr -Ofault.nc

# 转换结果可直接用于断层网格和逐行时间函数合成
grt fault2nc -Fconverted_param.inr -Ofault.nc
grt fault2nc -Fconverted_fsp.inr -Ofault.nc
grt lamb -H8/4.62/3.3 -N16/1 -Cconverted_param.inr+i0/0 -R200 -A0 -Dr0 -I1 -Osyn

# 输入格式必须为 p 或 f
expect_fail grt faultconvert -Fx/basic_inversion.param -C-8.3101/121.3517 -Oconverted_param.inr
# 未指定参考坐标
expect_fail grt faultconvert -Fp/basic_inversion.param -Oconverted_param.inr
# 参考纬度不能为极点
expect_fail grt faultconvert -Fp/basic_inversion.param -C90/0 -Oconverted_param.inr
# 段号须为正整数，且不可超过断层段数
expect_fail grt faultconvert -Fp/basic_inversion.param -C-8.3101/121.3517 -N0 -Oconverted_param.inr
expect_fail grt faultconvert -Ff/complete_inversion.fsp -C-8.3101/121.3517 -N2 -Oconverted_fsp.inr
# 输出不可覆盖输入
expect_fail grt faultconvert -Fp/basic_inversion.param -C-8.3101/121.3517 -Obasic_inversion.param
# 子断层数量须为正整数
sed 's/nx(Along-strike)= 25/nx(Along-strike)= -25/' basic_inversion.param > invalid.param
expect_fail grt faultconvert -Fp/invalid.param -C-8.3101/121.3517 -Oconverted_param.inr
sed 's/nx(Along-strike)= 25/nx(Along-strike)= 2.5/' basic_inversion.param > invalid.param
expect_fail grt faultconvert -Fp/invalid.param -C-8.3101/121.3517 -Oconverted_param.inr
# 段头尺寸必须带 km 单位
sed 's/Dy =  5.00km/Dy =  5.00/' basic_inversion.param > invalid.param
expect_fail grt faultconvert -Fp/invalid.param -C-8.3101/121.3517 -Oconverted_param.inr
# 元数据中的数值不可包含多余字符
sed 's/Dx = 10 km/Dx = 10abc km/' complete_inversion.fsp > invalid.fsp
expect_fail grt faultconvert -Ff/invalid.fsp -C-8.3101/121.3517 -Oconverted_fsp.inr
# 实际记录数必须覆盖完整网格
sed '$d' complete_inversion.fsp > invalid.fsp
expect_fail grt faultconvert -Ff/invalid.fsp -C-8.3101/121.3517 -Oconverted_fsp.inr
# FSP 总时长必须为正，不将零时长替换为脉冲
awk '!/^%/ && NF && !changed { $9 = 0; changed = 1 } { print }' complete_inversion.fsp > invalid.fsp
expect_fail grt faultconvert -Ff/invalid.fsp -C-8.3101/121.3517 -Oconverted_fsp.inr

rm -f basic_inversion.param complete_inversion.fsp legacy.param legacy.fsp legacy.inr invalid.param invalid.fsp count.fsp
rm -f converted_param.inr converted_fsp.inr reference_param.inr reference_fsp.inr fault.nc
rm -rf syn
