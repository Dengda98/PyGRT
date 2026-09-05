:author: 朱邓达
:date: 2026-09-01

.. include:: common_OPTs.rst_


lamb3
==============

:简介: 使用广义闭合解求解第三类 Lamb 问题

语法
-----------

**grt lamb3**
|-P|\ *nu*
|-T|\ *t1/t2/dt*
|-R|\ *dist*
|-D|\ *depsrc/deprcv*
|-A|\ *azimuth*
[ |-S|\ *+s<source-path>+r<receiver-path>* ]
[ **-h** ]


描述
--------

第三类 Lamb 问题指震源和观测点均位于半空间内部的情形。 
**lamb3** 模块实现的理论基础来源于《地震学中的 Lamb 问题（下）》第 8 章。
结果为与阶跃函数卷积后的无量纲位移，输出到标准输出，
第一列为无量纲时间，随后为 9 个 :math:`G_{ij}` 分量，按行优先顺序排列。
**-S** 指定的每个导数文件第一列同样为无量纲时间，随后为对应的
27 个导数分量，均按行优先顺序排列。所有输出均为无量纲量，实际物理量的恢复关系与
:doc:`lamb2` 模块中介绍的一致。

必选选项
-----------------

.. _-P:

**-P**\ *nu*
    半空间的泊松比 *nu*，要求范围在 (0, 0.5)。当 *nu* 距任一边界小于
    :math:`10^{-3}` 时给出数值稳定性警告，计算很可能失败

.. _-T:

**-T**\ *t1/t2/dt*
    无量纲时间序列 :math:`\bar{t}`，其中开始时间 *t1*、结束时间 *t2* 和时间间隔 *dt* 
    均以 :math:`\bar{t}` 为单位。
    :math:`\bar{t} = \dfrac{t}{T_S} = \dfrac{t}{r/\beta} = \dfrac{\beta t}{r}`，其中
    :math:`T_S = \dfrac{r}{\beta}` 是 S 波传播时间尺度，
    :math:`r` 为震源到接收点的直线距离，:math:`\beta` 为 S 波速度。

.. _-R:

**-R**\ *dist*
    源点到接收点的水平震中距。
    当 :math:`R/r'\leq 10^{-2}` 时给出数值稳定性警告，计算很可能失败，其中
    :math:`r'=\sqrt{R^2+(\mathit{depsrc}+\mathit{deprcv})^2}` 。

.. _-D:

**-D**\ *depsrc/deprcv*
    按 *depsrc/deprcv* 同时指定源点和观测点深度，二者都要求严格大于零。当任一
    深度与两点直线距离 :math:`r` 的比值小于 :math:`10^{-3}` 时，给出接近自由表面的
    数值稳定性警告，计算很可能失败。

.. _-A:

**-A**\ *azimuth*
    方位角，单位为度，要求范围在 [0, 360]

可选选项
--------------

.. _-S:

**-S**\ *+s<source-path>+r<receiver-path>*
    将源点坐标导数写入 *+s<source-path>* 指定的文件，将接收点坐标导数写入
    *+r<receiver-path>* 指定的文件。两个子选项可以只指定一个。每个文件包含
    无量纲时间列和 27 个对应导数列

参考文献
--------------

+ Feng, X., Zhang, H., 2021. Exact closed-form solutions for lamb’s problem—III: the case for buried source and receiver. Geophys. J. Int. 224, 517–532. https://doi.org/10.1093/gji/ggaa485
+ 张海明，冯禧，2024. 地震学中的Lamb问题（下）[M]. 北京：科学出版社.


示例
-------

详见 :doc:`/Lamb_problem/lamb3` 。
