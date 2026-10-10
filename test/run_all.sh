#!/bin/bash
# 默认运行全部测试，--local 仅本地，--remote 仅远程

set -euo pipefail

case "${1:-}" in
    ""|--local|--remote) ;;
    *) echo "用法: bash run_all.sh [--local|--remote]" >&2; exit 1 ;;
esac

# 下载远程正确性测试
if [[ "${1:-}" != --local ]]; then
    remote_repo=$(mktemp -d)
    remote_dirs=(_correct_full_wave _correct_static_disp _correct_surface_wave)
    git clone --depth 1 --filter=blob:none --sparse --branch main "https://github.com/Dengda98/dengda98.github.io.git" "$remote_repo"
    git -C "$remote_repo" sparse-checkout set assets/pygrt-tests
    for dir in "${remote_dirs[@]}"; do
        test -f "$remote_repo/assets/pygrt-tests/$dir/correctness.sh"
    done
    cp -R "$remote_repo/assets/pygrt-tests/." .
fi

if [[ "${1:-}" == --remote ]]; then
    test_dirs=("${remote_dirs[@]}")
else
    test_dirs=(*/)
fi

for dir in "${test_dirs[@]}"; do
    # 本地模式跳过已下载的远程测试目录
    [[ "${1:-}" == --local && "$dir" == _* ]] && continue
    (
        cd "$dir"
        for script in *.sh; do
            [[ -f "$script" ]] || continue
            echo "@ --------------- $dir, $script -------------"
            bash "$script"
        done
    )
done

# 删除下载的远程测试和临时检出目录
if [[ "${1:-}" != --local ]]; then
    rm -rf "$remote_repo" "${remote_dirs[@]}" compare_sac.py compare_nc.py
fi
