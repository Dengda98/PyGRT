#!/bin/bash

# 发起运行相关脚本

set -euo pipefail

rm -rf *.nc *.geo *.png *.tar.gz

head -n 10 coulomb-fault.inp > coulomb-fault.inp.head
echo "..." >> coulomb-fault.inp.head

bash run_grn.sh

bash run_syn.sh

bash run_gmt.sh

cp layer_disp.png cover.png

ex=$(basename $(pwd))
cd .. && tar -czvf ${ex}.tar.gz ${ex} && mv ${ex}.tar.gz ${ex} && cd -

rm -rf *.nc *.geo


