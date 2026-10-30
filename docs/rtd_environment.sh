#!/bin/bash
set -euo pipefail

uname -a
command -v python gcc gmt grt || true
python --version
gcc --version
getconf _NPROCESSORS_ONLN
if command -v lscpu >/dev/null; then
    lscpu
fi
for file in /sys/fs/cgroup/cpu.max /sys/fs/cgroup/cpu.stat /sys/fs/cgroup/memory.max /sys/fs/cgroup/cpu/cpu.cfs_quota_us /sys/fs/cgroup/cpu/cpu.cfs_period_us; do
    if [[ -f "$file" ]]; then
        echo "$file"
        cat "$file"
    fi
done
for name in OMP_NUM_THREADS OMP_WAIT_POLICY OPENBLAS_NUM_THREADS MKL_NUM_THREADS; do
    printf '%s=%s\n' "$name" "${!name-unset}"
done
python -m pip freeze
