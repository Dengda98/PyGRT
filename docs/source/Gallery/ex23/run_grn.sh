#!/bin/bash

set -euo pipefail

# USGS 配套 FSP 中的六层速度模型，地表接收点，仅计算剪切源及其空间导数
# 使用自适应 Filon 积分加快计算速度
rm -rf GRN
grt greenfn -Musgs.model -Ds2/35/3 -Dr0 -R180/300/10 -N500/1.0+a+n10 -L+a1e-2 -Gs -e -OGRN
