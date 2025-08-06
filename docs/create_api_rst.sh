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

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

.. toctree::
   :maxdepth: 2

EOF

for item in $(ls $SDIR); do
    echo $item
    if [ -d $SDIR/$item ]; then
        mkdir -p $PDIR/$item
        content=(
            "$item"
            "========================================="
            ""
            ".. toctree::"
            "   :maxdepth: 2"
            ""
        )
        printf "%s\n" "${content[@]}" > $PDIR/$item.rst

        for hnm in $(ls $SDIR/$item); do
            content=(
                "$hnm"
                "--------------------------------------"
                ""
                ".. doxygenfile:: $hnm"
                "   :project: h_PyGRT"
                ""
            )
            printf "%s\n" "${content[@]}" > $PDIR/$item/${hnm%%.*}.rst
            echo "   $item/${hnm%%.*}" >> $PDIR/$item.rst
        done

    elif [ -f $SDIR/$item ]; then
        content=(
            "$item"
            "--------------------------------------"
            ""
            ".. doxygenfile:: $item"
            "   :project: h_PyGRT"
            ""
        )
        printf "%s\n" "${content[@]}" > $PDIR/${item%%.*}.rst
    fi
    
    echo "   C_extension/${item%%.*}" >> $c_api_rst

    
done



# 根据源代码中pygrt/路径下的.py文件，自动在source/API/pygrt文件夹下生成.rst文件
PDIR="$APIDIR/pygrt"
rm -rf $PDIR
mkdir -p $PDIR 

SDIR="../pygrt"

py_api_rst="$APIDIR/py_api.rst"

cat > $py_api_rst <<EOF
Python API
========================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

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
