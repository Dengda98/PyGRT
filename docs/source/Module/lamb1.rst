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
[ |-C|\ *cbar* ]
[ |-Q|\ *phases* ]
[ **-h** ]


描述
--------

第一类 Lamb 问题指在半空间中，当源点和场点均位于地表时，求场点记录到的位移。
**lamb1** 模块的固定源和运动源部分分别参考了《地震学中的 Lamb 问题（下）》第 6 章和第 9 章的广义闭合解。
未指定 **-C** 时，固定源结果为与阶跃函数卷积后的无量纲位移，标准输出第一列为无量纲时间，
随后按行优先顺序排列 9 个 :math:`G_{ij}` 分量。指定 **-C<cbar>** 时启用垂直向下点力源的运动模式，
标准输出仅含无量纲时间和三个未与阶跃函数卷积的无量纲位移 :math:`\bar{u}_1`、:math:`\bar{u}_2`、
:math:`\bar{u}_3`。

为了得到实际物理单位的解，可进行：

.. math::

   G^H = \frac{\bar{G}^H}{\pi^2\mu r}

其中 :math:`r` 为震源到接收点的距离，:math:`\mu` 为剪切模量。运动源三个位移分量使用
相同归一化。

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

可选选项
-----------------

.. _-C:

**-C**\ *cbar*
    指定沿 :math:`x_1` 正方向匀速运动的垂直向下点力源无量纲速度 :math:`\bar{c}`，即
    :math:`c/\beta`，其中 :math:`c` 为实际速度，:math:`\beta` 为 S 波速度。
    要求 :math:`0 < c/\beta < v_R/\beta`，其中 :math:`v_R` 为 Rayleigh 波速度。
    此时台站不能位于 :math:`x_1` 轴上。此模式只输出时间和三个位移分量；

.. _-Q:

**-Q**\ *phases*
    固定源模式下仅保留指定震相项。*phases* 为由逗号分隔的 *P*、*S* 和 *R*。
    运动源模式不使用震相分解。无效震相名和重复震相名会给出警告，重复项只记录一次；
    剔除无效震相名后若没有有效震相，给出警告并输出全零波形。

参考文献
--------------

+ Feng, X., Zhang, H., 2018. Exact closed-form solutions for lamb’s problem. Geophys. J. Int. 214, 444–459. https://doi.org/10.1093/gji/ggy131
+ Feng, X., Zhang, H., 2020. Exact closed-form solutions for Lamb’s problem—II: a moving
  point load. *Geophys. J. Int.* 223, 1446–1459. https://doi.org/10.1093/gji/ggaa380
+ 张海明，冯禧，2024. 地震学中的Lamb问题（下）[M]. 北京：科学出版社.


示例
-------

详见 :doc:`/Lamb_problem/lamb1` 。
