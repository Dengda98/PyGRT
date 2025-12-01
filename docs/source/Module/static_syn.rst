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
[ |-F|\ *fn/fe/fz* ]
[ |-M|\ *strike/dip/rake* ]
[ |-T|\ *Mxx/Mxy/Mxz/Myy/Myz/Mzz* ]
[ |-N| ]
[ **-e** ]
[ **-h** ]


描述
--------

关于输出分量、震源机制的部分与 :doc:`syn` 模块相同，但输入和输出均使用 NetCDF 网格。


必选选项
----------

.. _-G:
    
**-G**\ *ingrid*
    存有静态格林函数的网格文件。

.. include:: explain_-S.rst_

.. include:: explain_-Ogrid.rst_



可选选项
--------

.. include:: explain_src.rst_

.. include:: explain_rot2ZNE.rst_

.. include:: explain_-esyn.rst_

.. include:: explain_-h.rst_



示例
-------

详见教程：

+ :doc:`/Tutorial/static/static_syn`