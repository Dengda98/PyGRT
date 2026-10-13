#!/bin/bash

set -euo pipefail

# 固定 ex19 使用的 USGS 产品版本，避免不同更新版本的结果混用
url="https://earthquake.usgs.gov/product/finite-fault/us6000tkt2_1/us/1787347333870"
for file in basic_inversion.param moment_rate.mr; do
    if [[ ! -f "$file" ]]; then
        curl -fLsS --retry 3 "$url/$file" -o "$file"
    fi
done

grt faultconvert -Fp/basic_inversion.param -C-8.3101/121.3517 -Ofault.inr
