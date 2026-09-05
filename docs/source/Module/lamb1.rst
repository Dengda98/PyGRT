:author: 朱邓达
:date: 2025-11-24

.. include:: common_OPTs.rst_


lamb1
==============

:简介: 使用广义闭合解求解第一类 Lamb 问题

语法
-----------

**grt lamb1**
|-P|\ *nu*
|-T|\ *t1/t2/dt*
|-A|\ *azimuth*
[ **-h** ]


描述
--------

第一类 Lamb 问题指在半空间中，当源点和场点均位于地表时，求场点记录到的位移。
**lamb1** 模块实现的理论基础来源于《地震学中的 Lamb 问题（下）》第 6 章。
结果为与阶跃函数卷积后的无量纲位移，输出到标准输出，
第一列为无量纲时间，随后为 9 个 :math:`G_{ij}` 分量，按行优先顺序排列。
实际的物理 Green 函数需要将输出结果除以 :math:`\pi^2\mu r`，即

.. math::

   G^H=\frac{\bar{G}^H}{\pi^2\mu r}

其中 :math:`r` 为震源到接收点的距离，:math:`\mu` 为剪切模量.

必选选项
-----------------

.. _-P:

**-P**\ *nu*
    半空间的泊松比 *nu*，要求范围在 (0, 0.5)。

.. _-T:

**-T**\ *t1/t2/dt*
    无量纲时间序列 :math:`\bar{t}`，其中开始时间 *t1*、结束时间 *t2* 和时间间隔 *dt* 
    均以 :math:`\bar{t}` 为单位。
    :math:`\bar{t} = \dfrac{t}{T_S} = \dfrac{t}{r/\beta} = \dfrac{\beta t}{r}`，其中
    :math:`T_S = \dfrac{r}{\beta}` 是 S 波传播时间尺度，
    :math:`r` 为震源到接收点的直线距离，:math:`\beta` 为 S 波速度。

.. _-A:

**-A**\ *azimuth*
    方位角，单位为度。

参考文献
--------------

+ Feng, X., Zhang, H., 2018. Exact closed-form solutions for lamb’s problem. Geophys. J. Int. 214, 444–459. https://doi.org/10.1093/gji/ggy131
+ 张海明，冯禧，2024. 地震学中的Lamb问题（下）[M]. 北京：科学出版社.


示例
-------

详见 :doc:`/Lamb_problem/lamb1` 。
