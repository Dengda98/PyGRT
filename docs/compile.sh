#!/bin/bash

set -euo pipefail

chmod +x *.sh

# 先在线构建PyGRT
cd ../pygrt/C_extension 
make
cd -

# 添加到环境变量
# cp ../pygrt/C_extension/bin/* /usr/local/bin/
# export PATH=$(realpath ../pygrt/C_extension/bin):$PATH
# echo "export PATH=$(realpath ../pygrt/C_extension/bin):$PATH" >> ~/.bashrc
cp ../pygrt/C_extension/bin/* ${READTHEDOCS_VIRTUALENV_PATH}/bin/
echo "-------------------------"
echo $PATH
# echo "-------------------------"
# echo $(ls /usr/local/bin/* -l)
echo "-------------------------"
echo ${READTHEDOCS_VIRTUALENV_PATH}
echo ${READTHEDOCS_REPOSITORY_PATH}
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

# 构建
if [[ $1 == '1' ]]; then
make html
fi
# sphinx-autobuild -j auto --port 8000 source build