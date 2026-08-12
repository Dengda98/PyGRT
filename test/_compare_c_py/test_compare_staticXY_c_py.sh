#!/bin/bash
# 比较静态合成指定新 XY 网格时 C CLI 与 Python API 的结果

set -euo pipefail

cd "$(dirname "$0")"

ROOT="$(cd ../.. && pwd)"
export PATH="${ROOT}/pygrt/C_extension/bin:${PATH}"

echo "grt = $(command -v grt)"
echo "python = $(command -v python)"

echo "========== staticXY C vs Python =========="
python -u compare_staticXY.py

# 清理测试产生的临时文件，保持目录整洁
rm -rf _work_compare_staticXY __pycache__

echo "All tests in test_compare_staticXY_c_py.sh passed."
