#!/bin/bash

# 根据源代码中include/路径下的.h文件，自动在source/API/C_extension文件夹下生成.rst文件

APIDIR="source/API"

PDIR="$APIDIR/C_extension"
rm -rf $PDIR
mkdir -p $PDIR 

SDIR="../pygrt/C_extension/include"

c_api_rst="$APIDIR/c_api.rst"

cat > $c_api_rst <<EOF
C API
===============================

.. toctree::
   :maxdepth: 2

EOF

# 递归处理目录和文件
process_directory() {
    local source_dir=$1
    local target_dir=$2
    local relative_path=$3

    echo $source_dir $target_dir $relative_path
    
    for item in $(ls $source_dir); do
        local source_item="$source_dir/$item"
        local target_item="$target_dir/$item"
        local relative_item="$relative_path/$item"
        
        if [ -d $source_item ]; then
            # 创建对应的目标目录
            mkdir -p "$target_item"
            
            # 为目录创建.rst文件
            local rst_file="${target_dir}/${item}.rst"
            content=(
                "$item"
                "========================================="
                ""
                ".. toctree::"
                "   :maxdepth: 2"
                ""
            )
            printf "%s\n" "${content[@]}" > "$rst_file"

            # 将目录添加到父目录的toctree中
            if [ "$relative_path" != "." ]; then
                local parent_rst="${target_dir}.rst"
                echo "   $(basename $source_dir)/${item}" >> "$parent_rst"
            else
                # 如果是顶级目录下的目录，添加到主toctree
                echo "   C_extension/${item}" >> "$c_api_rst"
            fi
            
            # 递归处理子目录
            process_directory "$source_item" "$target_item" "$relative_item"
            
        elif [ -f $source_item ] && [[ $item == *.h ]]; then
            # 处理.h文件
            # local base_name="${item%%.*}"
            local base_name="${item}"
            content=(
                "$item"
                "--------------------------------------"
                ""
                ".. doxygenfile:: ${relative_item#./}"
                "   :project: h_PyGRT"
                ""
            )
            printf "%s\n" "${content[@]}" > "${target_dir}/${base_name}.rst"
            
            # 如果是顶级目录下的文件，添加到主toctree
            if [ "$relative_path" == "." ]; then
                echo "   C_extension/${base_name}" >> "$c_api_rst"
            else
                # 将文件添加到父目录的toctree中
                local parent_rst="${target_dir}.rst"
                echo "   $(basename ${target_dir})/${base_name}" >> "$parent_rst"
            fi
        fi
    done
}

# 从SDIR开始处理
process_directory "$SDIR" "$PDIR" "."




# 根据源代码中pygrt/路径下的.py文件，自动在source/API/pygrt文件夹下生成.rst文件
PDIR="$APIDIR/pygrt"
rm -rf $PDIR
mkdir -p $PDIR 

SDIR="../pygrt"

py_api_rst="$APIDIR/py_api.rst"

cat > $py_api_rst <<EOF
Python API
========================

.. toctree::
   :maxdepth: 2

EOF

for py in $(ls -1 $SDIR/*.py | xargs -n1 basename); do
    if [[ "$py" == "__init__.py" || "$py" == "_version.py" || "$py" == "print.py" ]]; then
        continue
    fi

    content=(
        "pygrt.${py%%.*}"
        "--------------------------------------"
        ""
        ".. automodule:: pygrt.${py%%.*}"
        "   :members:"
        "   :undoc-members:"
        "   :show-inheritance:"
        ""
    )

    printf "%s\n" "${content[@]}" > $PDIR/${py%%.*}.rst
    echo "   pygrt/${py%%.*}" >> $py_api_rst

done
