#!/bin/bash

set -euo pipefail

# 计算格林函数库
grt static greenfn -Mprem.flat.20 -Ds1/36/2 -Dr0 -R0/500/2.0 -Ostgrn.nc -e

