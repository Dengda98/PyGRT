#!/bin/bash

# 根据源代码中include/路径下的.h文件，自动在source/API/C_extension文件夹下生成.rst文件

APIDIR="source/API"

PDIR="$APIDIR/C_extension"
rm -rf $PDIR
mkdir -p $PDIR 

SDIR="../pygrt/C_extension/include"

c_api_rst="$APIDIR/c_api.rst"

cat > $c_api_rst <<EOF
C extension API
===============================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

.. toctree::
   :maxdepth: 2

EOF

for dir1 in $(ls $SDIR); do
    mkdir -p $PDIR/$dir1
    cat > $PDIR/$dir1.rst <<EOF
$dir1
=========================================

.. toctree::
   :maxdepth: 2

EOF
    for hnm in $(ls $SDIR/$dir1); do
        cat > $PDIR/$dir1/${hnm%%.*}.rst <<EOF
$hnm
--------------------------------------

.. doxygenfile:: $hnm
    :project: h_PyGRT
EOF

    echo "   $dir1/${hnm%%.*}" >> $PDIR/$dir1.rst
    done

    echo "   C_extension/$dir1" >> $c_api_rst
done



# 根据源代码中pygrt/路径下的.py文件，自动在source/API/pygrt文件夹下生成.rst文件
PDIR="$APIDIR/pygrt"
rm -rf $PDIR
mkdir -p $PDIR 

SDIR="../pygrt"

py_api_rst="$APIDIR/py_api.rst"

cat > $py_api_rst <<EOF
**PyGRT** package API
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

    cat > $PDIR/${py%%.*}.rst <<EOF
pygrt.${py%%.*}
--------------------------------------

.. automodule:: pygrt.${py%%.*}
    :members:
    :undoc-members:
    :show-inheritance:

EOF

    echo "   pygrt/${py%%.*}" >> $py_api_rst

done
