#!/bin/bash

set -euo pipefail

# 使用 USGS 配套 FSP 中的六层模型计算静态格林函数库
grt static greenfn -Musgs.model -Ds1/35/2 -Dr0 -R0/500/2.0 -Ostgrn.nc -e

