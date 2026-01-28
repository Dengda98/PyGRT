:author: 朱邓达
:date: 2025-09-22

.. include:: common_OPTs.rst_


static_syn
==================

:简介: 指定震源机制，根据静态格林函数合成三分量位移（及其空间导数）


语法
-----------

**grt static syn** 
|-G|\ *ingrid*
|-S|\ [**u**]\ *scale*
|-O|\ *outgrid*
[ |-X|\ *x1/x2/dx* ]
[ |-Y|\ *y1/y2/dy* ]
[ |-F|\ *fn/fe/fz* ]
[ |-M|\ *strike/dip/rake* ]
[ |-T|\ *Mxx/Mxy/Mxz/Myy/Myz/Mzz* ]
[ |-N| ]
[ **-e** ]
[ **-h** ]


描述
--------

关于输出分量、震源机制的部分与 :doc:`syn` 模块相同，但输入和输出均使用 NetCDF 网格。

**static syn** 模块计算时使用的 XY 网格默认继承 :doc:`static_greenfn` 模块设置的 **-X** 和 **-Y** 。
也支持设置 |-X| 和 |-Y| 来指定新的 XY 网格，此时每个节点将使用最近震中距的格林函数来进行合成，
因此这是一种近似方法，程序会打印震中距之差的统计信息，以方便评估近似效果。
但近似好处是一旦静态格林函数计算好，合成阶段使用新的网格也可以复用格林函数，无需再计算。


必选选项
----------

.. _-G:
    
**-G**\ *ingrid*
    存有静态格林函数的网格文件。

.. include:: explain_-S.rst_

.. include:: explain_-Ogrid.rst_



可选选项
--------

.. include:: explain_-XYgrid.rst_

.. include:: explain_src.rst_

.. include:: explain_rot2ZNE.rst_

.. include:: explain_-esyn.rst_

.. include:: explain_-h.rst_



示例
-------

详见教程：

+ :doc:`/Tutorial/static/static_syn`