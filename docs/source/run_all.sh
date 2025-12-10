#!/bin/bash

set -euo pipefail

# 运行教程中所有示例的脚本，方便将成图过程放在readthedocs服务器上
# 为方便统一发起任务，本脚本将在本目录下寻找所有名字为 run.sh 的脚本，
# 切换到对应目录下，执行 bash run.sh


fname="run.sh"
for path in $(find -type f -name "$fname"); do
    cd $(dirname $path)
    chmod +x $fname
    echo "--------------------------- $fname -----------------------------"
    bash $fname 
    cd -
done
