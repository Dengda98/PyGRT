:author: 朱邓达
:date: 2026-08-24

.. include:: common_OPTs.rst_


static_coulomb
========================

:简介: 根据投影后的法向应力和剪应力计算库伦应力变化


语法
-----------

**grt static coulomb**
|-G|\ *ingrid*
|-F|\ *friction*
[ **-h** ]


描述
--------

**static coulomb** 模块读取输入静态合成文件中的 **sigma_n** 和 **tau_s**，
按照下式计算每个接收点上的库伦应力变化：

.. math::

    \Delta CFS = \Delta \tau_\text{s} + \mu^{'} \times \Delta \sigma_\text{n}

其中 **friction** 为无量纲等效摩擦系数。**sigma_n** 和 **tau_s** 由
:doc:`static_sproj` 模块写入。结果变量 **coulomb** 与输入应力使用相同的单位，
即 dyne/cm²（0.1 Pa）。


必选选项
----------

.. _-G:

**-G**\ *ingrid*
    输入静态合成文件。文件必须包含 **sigma_n** 和 **tau_s** 变量，
    应先执行 :doc:`static_sproj` 模块。

.. _-F:

**-F**\ *friction*
    无量纲等效摩擦系数，必须为有限的非负数。


可选选项
----------

.. include:: explain_-h.rst_
