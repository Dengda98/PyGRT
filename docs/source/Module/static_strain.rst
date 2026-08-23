:author: 朱邓达
:date: 2025-09-22

.. include:: common_OPTs.rst_


static_strain
=================

:简介: 根据合成的静态位移空间导数计算静态应变张量


语法
-----------

**grt static strain** *ingrid* [ **-h** ]

描述
--------------

**static strain** 模块支持读取 :doc:`static_syn` 的输出文件，计算静态应变张量，并将结果写回原文件。
合成时需使用 **-e**。公式为

.. math:: 

    e_{ij} = \dfrac{1}{2} \left( u_{i,j} + u_{j,i} \right) = \dfrac{1}{2} \left( \dfrac{\partial u_i}{\partial x_j} + \dfrac{\partial u_j}{\partial x_i}  \right)

参数 *ingrid* 为 :doc:`static_syn` 的输出文件。

**static strain** 模块将六个分量写入 *ingrid* 文件，变量名为 ``strain_??``，
其中 ``??`` 代表六个分量名，即上述公式中的下标 :math:`ij` 。
如果合成的位移使用 ZRT 分量，则六个分量分别为 *ZZ,ZR,ZT,RR,RT,TT* ；
如果合成的位移使用 ZNE 分量，则六个分量分别为 *ZZ,ZN,ZE,NN,NE,EE* 。



示例
-------

详见教程：

+ :doc:`/Tutorial/dynamic/strain_stress`
