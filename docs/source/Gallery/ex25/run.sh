#!/bin/bash

set -euo pipefail

python 01_create_input_files.py

# 本例计算量较大，文档构建直接使用随示例保存的图件
# 如需重新计算，在本目录下手动执行 bash reproduce.sh

cp ambient_figures/profiles/STA001.png cover.png

ex=$(basename "$(pwd)")
rm -f "$ex.tar.gz"
cd .. && tar -czvf "${ex}.tar.gz" "${ex}" && mv "${ex}.tar.gz" "${ex}" && cd -
