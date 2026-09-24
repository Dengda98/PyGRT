:author: 朱邓达
:date: 2025-11-14

第一类 Lamb 问题（lamb1 模块）
==================================

第一类 Lamb 问题是指，在半空间模型中，源点和场点均位于地表，求解场点记录到的位移。

.. tabs::  

   .. group-tab:: CLI
      
      :command:`grt` 命令提供了模块 :doc:`/Module/lamb1` 求解第一类 Lamb 问题。

      .. literalinclude:: run/run.sh
         :language: bash
         :start-after: BEGIN LAMB1
         :end-before: END LAMB1

      使用重定向将结果保存到文件 *lamb1.txt* 中，其内容格式如下，
      记录了无量纲时间 :math:`\bar{t}` 和无量纲 Green 函数 :math:`G_{ij}` ，
      表示 :math:`j` 方向的力激发的 :math:`i` 方向的位移。

      .. literalinclude:: run/head_lamb1
         :language: text

   .. group-tab:: Python    

      Python 提供了函数 :func:`lamb1() <pygrt.utils.lamb1>` 求解第一类 Lamb 问题。
      函数返回无量纲 Green 函数数组 ``G``

      .. literalinclude:: run/lamb1_plot_time.py
         :language: python
         :start-after: BEGIN LAMB1
         :end-before: END LAMB1

书中结果复现
--------------------

:download:`lamb1_plot_time.py <run/lamb1_plot_time.py>`

.. figure:: run/lamb1_time.svg
   :align: center

   复现了原书中的图 6.6.4


频域解和时域解的对比
-------------------------------

对比观察可发现频域解在波形突变出有明显的 Gibbs 效应。

:download:`lamb1_plot_freq_time.py <run/lamb1_plot_freq_time.py>`

.. figure:: run/lamb1_compare_freq_time.svg
   :align: center

运动源
--------------------

基于《地震学中的 Lamb 问题（下）》第 9 章给出的广义闭合解，
可以计算出沿 :math:`x_1` 正方向匀速移动的垂直向下的点力源所激发的三分量位移，
速度以 :math:`c/\beta` 指定。

.. tabs::

   .. group-tab:: CLI

      :doc:`/Module/lamb1` 通过 **-C<cbar>** 来指定垂直力源的无量纲移动速度。

      .. literalinclude:: run/run.sh
         :language: bash
         :start-after: BEGIN LAMB1 MOVING
         :end-before: END LAMB1 MOVING

      使用重定向将结果保存到文件 *lamb1_moving.txt* 中，其内容格式如下，
      记录了无量纲时间 :math:`\bar{t}` 和无量纲三分量位移。

      .. literalinclude:: run/head_lamb1_moving
         :language: text

   .. group-tab:: Python

      Python 中可在 :func:`lamb1() <pygrt.utils.lamb1>` 函数中指定 ``cbar``
      参数来指定垂直力源的无量纲移动速度。

      .. literalinclude:: run/lamb1_plot_time.py
         :language: python
         :start-after: BEGIN LAMB1 MOVING
         :end-before: END LAMB1 MOVING

以下复现书中第九章的一些结果。

:download:`lamb1_plot_moving_source_9_4_1.py <run/lamb1_plot_moving_source_9_4_1.py>`

.. figure:: run/lamb1_moving_source_9_4_1.svg
   :align: center

   复现图 9.4.1：固定垂直点力与低速运动垂直点力的位移分量比较。

:download:`lamb1_plot_moving_source_9_4_4.py <run/lamb1_plot_moving_source_9_4_4.py>`

.. figure:: run/lamb1_moving_source_9_4_4.svg
   :align: center

   复现图 9.4.4：沿 :math:`x_1` 轴两侧测线的地表位移波形。

:download:`lamb1_plot_moving_source_9_4_5.py <run/lamb1_plot_moving_source_9_4_5.py>`

.. figure:: run/lamb1_moving_source_9_4_5.svg
   :align: center

   复现图 9.4.5：运动方向前方与后方测线的地表位移波形。
