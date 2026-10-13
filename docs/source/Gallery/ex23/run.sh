#!/bin/bash

set -euo pipefail

bash download.sh
# 计算仅运行 CLI 脚本，Python 版本用于参考展示
bash run_grn.sh
bash run_syn.sh
bash run_gmt.sh
python plot.py
python plot_moment_rate.py

cp slip_delay.svg cover.svg

ex=$(basename "$(pwd)")
rm -f "$ex.tar.gz"
cd .. && tar -czvf "${ex}.tar.gz" "${ex}" && mv "${ex}.tar.gz" "${ex}" && cd -
