:author: 朱邓达
:date: 2026-05-20

.. include:: common_OPTs.rst_


disp2asc
=================

:简介: 将 |NetCDF| 格式的面波频散结果转为 ASCII 文件


语法
-----------

**grt disp2asc** 
|-C|\ *path*
|-U|\ *path*
[ |-F|\ *f1*\ [/\ *f2*][/\ *df*][**+p**] ]
[ |-N|\ [*n1*][/\ *n2*][/\ *dn*] ]
[ **-h** ]

描述
--------

**disp2asc** 模块可将 :doc:`eigenv` 模块计算出的相速度频散结果转为方便阅读的 ASCII 文件，
结果输出到标准输出，包含四列::
    
    频率(Hz)      阶数      相速度(km/s)       久期函数层位

也可转换 :doc:`eigenfn` **-U** 模块计算出的群速度频散结果，结果包含四列::
    
    频率(Hz)      阶数      相速度(km/s)       群速度(km/s)



必选选项
----------

.. include:: explain_-Ceigv.rst_

.. _-U:

**-U**\ **path**
    :doc:`eigenfn` 模块使用 **-U** 计算的群速度频散结果保存文件路径。

.. warning::

    |-C| 和 |-U| 必须设置也仅能设置其中一个。


可选选项
--------

.. include:: explain_-Ffreqs.rst_

.. include:: explain_-Nmodes.rst_

.. include:: explain_-h.rst_



示例
-------
