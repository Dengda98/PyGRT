:author: 朱邓达
:date: 2025-09-22

.. include:: common_OPTs.rst_


static_strain
=================

:简介: 根据合成的静态位移空间导数计算静态应变张量


语法
-----------

**grt static strain** *syngrid*

描述
--------------

**static_strain** 模块计算静态应变张量，要求 :doc:`static_syn` 模块计算时使用 **-e** 以合成位移空间导数。公式为

.. math:: 

    e_{ij} = \dfrac{1}{2} \left( u_{i,j} + u_{j,i} \right) = \dfrac{1}{2} \left( \dfrac{\partial u_i}{\partial x_j} + \dfrac{\partial u_j}{\partial x_i}  \right)

参数 *syngrid* 表示 :doc:`static_syn` 模块中使用 **-O** 指定的输出网格文件。
**static_strain** 模块将合成的六个分量写入 *syngrid* 网格文件，变量名为 ``strain_??`` ,
其中 ``??`` 代表六个分量名，即上述公式中的下标 :math:`ij` 。
如果合成的位移使用 ZRT 分量，则六个分量分别为 *ZZ,ZR,ZT,RR,RT,TT* ；
如果合成的位移使用 ZNE 分量，则六个分量分别为 *ZZ,ZN,ZE,NN,NE,EE* 。



示例
-------

详见教程：

+ :doc:`/Tutorial/dynamic/strain_stress`
