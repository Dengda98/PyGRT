#!/bin/bash
# Run all tests

set -euo pipefail

# 下载其它脚本
remote_repo=$(mktemp -d)
git clone --depth 1 --filter=blob:none --sparse --branch main "https://github.com/Dengda98/dengda98.github.io.git" "$remote_repo"
git -C "$remote_repo" sparse-checkout set assets/pygrt-tests
remote_dirs=(_correct_full_wave _correct_static_disp _correct_surface_wave)
for dir in "${remote_dirs[@]}"; do
    test -f "$remote_repo/assets/pygrt-tests/$dir/correctness.sh"
done
cp -R "$remote_repo/assets/pygrt-tests/." .
rm -rf "$remote_repo"

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
rm -rf "${remote_dirs[@]}" compare_sac.py compare_nc.py
