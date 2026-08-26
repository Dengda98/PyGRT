#!/bin/bash

set -euo pipefail

chmod +x *.sh

# 先在线构建PyGRT
cd ../pygrt/C_extension 
make clean && make -j
cd -

# 将编译出的 grt 放入当前构建环境的 bin 目录
grt_source_dir="$(realpath ../pygrt/C_extension/bin)"
if [[ -n "${READTHEDOCS_VIRTUALENV_PATH:-}" ]]; then
    grt_bin_dir="${READTHEDOCS_VIRTUALENV_PATH}/bin"
elif [[ -n "${CONDA_PREFIX:-}" ]]; then
    grt_bin_dir="${CONDA_PREFIX}/bin"
else
    grt_bin_dir="${grt_source_dir}"
fi

if [[ "${grt_bin_dir}" != "${grt_source_dir}" ]]; then
    cp "${grt_source_dir}"/* "${grt_bin_dir}/"
fi
export PATH="${grt_bin_dir}:${PATH}"
echo "-------------------------"
echo "${PATH}"
# echo "-------------------------"
# echo $(ls /usr/local/bin/* -l)
echo "-------------------------"
echo "${grt_bin_dir}"
echo "${READTHEDOCS_REPOSITORY_PATH:-}"
# if [[ $(which grt) == "" ]]; then
# echo "export PATH=$(realpath ../pygrt/C_extension/bin):\$PATH" >> ~/.bashrc
# source ~/.bashrc
# fi
grt -h
# 使用PyGRT运行文档需要的示例结果
cd source && chmod +x *.sh && ./run_all.sh && cd -

# 清空构建的旧文档
make clean

# 生成api对应的.rst文件
./create_api_rst.sh

# 执行doxygen
doxygen doxyfile_h

# sphinx-autobuild -j auto --port 8000 source build
