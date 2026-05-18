#!/bin/bash
# Run all tests

set -euo pipefail

# 下载其它脚本
wget "https://raw.githubusercontent.com/Dengda98/dengda98.github.io/refs/heads/main/assets/pygrt-tests/pygrt-tests.tar.gz"
tar -xzvf pygrt-tests.tar.gz

dirs=$(find . -maxdepth 1 -mindepth 1 -type d | sort)

for dir in $dirs; do
    cd $dir > /dev/null
    for fname in $(ls *.sh); do
        echo "@ --------------- $dir, $fname -------------"
        bash $fname
    done
    cd - > /dev/null
done

# 删除下载的脚本
tar -tf pygrt-tests.tar.gz | awk -F/ '{print $1}' | uniq | xargs rm -rf
rm pygrt-tests.tar.gz