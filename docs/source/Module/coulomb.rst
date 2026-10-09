:author: 朱邓达
:date: 2026-10-03

.. include:: common_OPTs.rst_


coulomb
==========

:简介: 根据投影后的法向应力和剪应力计算动态库伦应力变化


语法
-----------

**grt coulomb**
|-G|\ *syn_dir*
|-F|\ *friction*
[ **-h** ]


描述
--------

:doc:`coulomb` 模块读取 :doc:`sproj` 生成的 **sigma_n.sac** 和 **tau_s.sac**，
按下式计算动态库伦应力变化：

.. math::

    \Delta CFS(t) = \Delta \tau_\text{s}(t) + \mu^{'} \times \Delta \sigma_\text{n}(t)

其中 **friction** 为无量纲等效摩擦系数。结果保存为输入目录中的 **coulomb.sac**，
单位与输入应力相同，即 dyne/cm²（0.1 Pa）。


必选选项
----------

.. _-G:

**-G**\ *syn_dir*
    包含 **sigma_n.sac** 和 **tau_s.sac** 的单台目录或多台根目录。
    多台结果逐接收目录计算并保存。

.. _-F:

**-F**\ *friction*
    无量纲等效摩擦系数，须为非负数。
