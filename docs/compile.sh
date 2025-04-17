#!/bin/bash

chmod +x *.sh

# 清空
make clean

# 生成api对应的.rst文件
./create_api_rst.sh

# 执行doxygen
doxygen doxyfile_h

# 构建
if [[ $1 == '1' ]]; then
make html
fi
# sphinx-autobuild --port 8000 source build