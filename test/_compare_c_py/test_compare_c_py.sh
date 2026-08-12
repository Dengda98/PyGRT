#!/bin/bash
# 比较 C CLI 与 Python 文件工作流：参数映射 + 动态/静态结果

set -euo pipefail

cd "$(dirname "$0")"

# 优先使用源码树内的 grt
ROOT="$(cd ../.. && pwd)"
export PATH="${ROOT}/pygrt/C_extension/bin:${PATH}"

echo "grt = $(command -v grt)"
echo "python = $(command -v python)"

echo "========== CLI argument mapping =========="
python -u test_cli_args.py

echo "========== end-to-end C vs Python =========="
python -u compare.py

# 清理测试产生的临时文件，保持目录整洁
rm -rf _work_compare _tmp_args_* _tmp_tensor_* __pycache__

echo "All tests in test_compare_c_py.sh passed."
