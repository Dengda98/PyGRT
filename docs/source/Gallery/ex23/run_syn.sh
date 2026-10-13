#!/bin/bash

set -euo pipefail

# 200 km 的地表接收点，不再额外细分子断层
# -e 同时输出九个位移空间导数
grt syn -GGRN -Cfault.inr+i0/0 -R200 -A10 -Dr0 -I1 -N -e -Osyn

# 法向应力张开为正，剪应力沿 rake 方向为正
grt stress  -Gsyn
grt sproj   -Gsyn -M90/32/90
grt coulomb -Gsyn -F0.4
