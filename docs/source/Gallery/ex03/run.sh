#!/bin/bash 

set -euo pipefail

rm -rf GRN* milrow_sdep*

# 程序受限，这里暂且上传静态图片
# bash run_milrow_grt.sh_
# bash run_milrow_cps.sh_
# python plot_cps_grt.py
# python plot_cps_pygrt.py

cp compare_cps_grt.svg cover.svg

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf GRN* milrow_sdep*