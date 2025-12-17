:author: 朱邓达
:date: 2025-12-10

.. include:: common_OPTs.rst_


kernel
=================

:简介: 使用广义反射透射系数矩阵法计算不同频率、不同相速度的核函数

语法
----------

**grt kernel**
|-M|\ *model*
|-D|\ *depsrc/deprcv*
|-F|\ *f1/f2/df*\ [**+w**\ *zeta*]
|-C|\ [*cmin/cmax/*]\ *dc*
|-O|\ *outdir*
[ |-P|\ *nthreads* ]
[ **-e** ]
[ **-h** ]


描述
--------

**kernel** 模块算是 :doc:`greenfn` 模块中 **-S** 选项的便捷版本。
通过直接指定离散的频率和相速度，计算对应的核函数。
与 :doc:`greenfn` **-S** 的输出相比， **kernel** 模块输出的是理论核函数，输出文件将以 ``C`` 开头，
表示采样点位置记录的是相速度。
输出文件为二进制格式，可使用 :doc:`ker2asc` 转为文本格式。


必选选项
----------

.. include:: explain_-M.rst_

.. include:: explain_-D.rst_

.. _-F:

**-F**\ *f1/f2/df*\ [**+w**\ *zeta*]
    频率范围 （Hz）。

    + *f1* - 最小频率，如果为 0 则自动跳到下一个 *df* 。
    + *f2* - 最大频率。
    + *df* - 频率间隔。

    **+w**\ *zeta* 用于在频率上增加一个微小虚部，即 :math:`\omega \leftarrow \omega - i \omega_I` ，
    其中 :math:`\omega_I = \zeta \cdot \pi \cdot \Delta f` 。 
    这里虚频率的效果是使核函数在极点附近的峰值范围变宽，方便可视化。[默认 *zeta=20*]

.. _-C:

**-C**\ [*cmin/cmax/*]\ *dc*
    相速度范围 （km/s）。

    + *cmin* - 最小相速度。
    + *cmax* - 最大相速度。
    + *dc*   - 相速度间隔。
    
    如果仅指定  **-C**\ *dc* ，则自动设置 *cmin* 为 0.8 倍的模型最小非零速度，
    *cmax* 为模型最大速度。

.. include:: explain_-O.rst_


可选选项
--------

.. include:: explain_-P.rst_

**-e**
    改为输出位移对z偏导的核函数。

.. include:: explain_-h.rst_


示例
-------

详见 :doc:`/Advanced/kernel/kernel` 。