:author: 朱邓达
:date: 2025-09-22

.. include:: common_OPTs.rst_


static_rotation
=======================

:简介: 根据合成的静态位移空间导数计算静态旋转张量


语法
-----------

**grt static rotation** *syngrid*

描述
--------------

**static_rotation** 模块计算静态旋转张量，要求 :doc:`static_syn` 模块计算时使用 **-e** 以合成位移空间导数。公式为

.. math:: 

    w_{ij} = \dfrac{1}{2} \left( u_{i,j} - u_{j,i} \right) = \dfrac{1}{2} \left( \dfrac{\partial u_i}{\partial x_j} - \dfrac{\partial u_j}{\partial x_i}  \right)

参数 *syngrid* 表示 :doc:`static_syn` 模块中使用 **-O** 指定的输出目录。
**static_rotation** 模块将合成的三个分量写入 *syngrid* 网格文件，变量名为 ``rotation_??`` ,
其中 ``??`` 代表三个分量名，即上述公式中的下标 :math:`ij` 。
如果合成的位移使用 ZRT 分量，则三个分量分别为 *ZR,ZT,RT* ；
如果合成的位移使用 ZNE 分量，则三个分量分别为 *ZN,ZE,NE* 。



示例
-------

详见教程：

+ :doc:`/Tutorial/dynamic/strain_stress`
