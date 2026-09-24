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
[ |-Q|\ *phases* ]
[ |-S|\ [*+s<source-path>*][*+r<receiver-path>*][*+m<mixed-path>*] ]
[ **-h** ]


描述
--------

第三类 Lamb 问题指震源和观测点均位于半空间内部的情形。 
**lamb3** 模块实现的理论基础来源于《地震学中的 Lamb 问题（下）》第 8 章。
结果为与阶跃函数卷积后的无量纲位移，输出到标准输出。标准输出第一列为无量纲时间，
随后按行优先顺序排列 9 个 :math:`G_{ij}` 分量。**-S** 指定的每个导数文件第一列同样为
无量纲时间，随后按行优先顺序排列 27 个一阶导数或 81 个二阶混合导数分量。
源点一阶导数的排列顺序为 :math:`k',i,j`，接收点一阶导数为 :math:`k,i,j`，
二阶混合导数的排列顺序为 :math:`k,k',i,j`。
二阶混合导数在闭合解内部由对应的时间积分项连续求两次时间导数得到，输出已经是最终结果。

为了得到实际物理单位的解，可进行：

.. math::

   G^H=\frac{\bar{G}^H}{\pi^2\mu r},\qquad
   G^H_{,k'}=\frac{\bar{G}^H_{,k'}}{\pi^2\mu r^2},\qquad
   G^H_{,k}=\frac{\bar{G}^H_{,k}}{\pi^2\mu r^2},\qquad
   G^H_{,k,k'}=\frac{\bar{G}^H_{,k,k'}}{\pi^2\mu r^3}

其中 :math:`r` 为震源到接收点的距离，:math:`\mu` 为剪切模量。

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
    深度与两点直线距离 :math:`r` 的比值小于等于 :math:`10^{-3}` 时，给出接近自由表面的
    数值稳定性警告，计算很可能失败。

.. _-A:

**-A**\ *azimuth*
    方位角，单位为度，要求范围在 [0, 360]

可选选项
--------------

.. _-Q:

**-Q**\ *phases*
    仅保留指定震相项。*phases* 为由逗号分隔的 *P*、*S*、*PP*、*SS*、*PS*、
    *SP* 和 *sPs*。其中 *P* 和 *S* 分别为直达
    P、S 波，*PP*、*SS*、*PS* 和 *SP* 为自由表面反射或转换波，
    *sPs* 为滑行波。

.. _-S:

**-S**\ [*+s<source-path>*][*+r<receiver-path>*][*+m<mixed-path>*]
    指定相关物理量的保存路径：
    至少指定一个子选项，多个子选项可以同时使用

    + *+s<source-path>* - 源点坐标一阶导数
    + *+r<receiver-path>* - 接收点坐标一阶导数
    + *+m<mixed-path>* - 接收点和源点坐标的二阶混合导数

参考文献
--------------

+ Feng, X., Zhang, H., 2021. Exact closed-form solutions for lamb’s problem—III: the case for buried source and receiver. Geophys. J. Int. 224, 517–532. https://doi.org/10.1093/gji/ggaa485
+ 张海明，冯禧，2024. 地震学中的Lamb问题（下）[M]. 北京：科学出版社.


示例
-------

详见 :doc:`/Lamb_problem/lamb3` 。
