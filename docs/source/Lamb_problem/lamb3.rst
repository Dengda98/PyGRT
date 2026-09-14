:author: 朱邓达
:date: 2026-09-01

第三类 Lamb 问题
===================

第三类 Lamb 问题中，源点和观测点都位于半空间内部。

.. tabs::

   .. group-tab:: CLI

      :command:`grt` 命令提供了模块 :doc:`/Module/lamb3` 求解第三类 Lamb 问题。

      .. literalinclude:: run_lamb3/run.sh
         :language: bash
         :start-after: BEGIN LAMB3
         :end-before: END LAMB3

      使用重定向将结果保存到文件 *lamb3.txt* 中，其内容格式如下，
      记录了无量纲时间 :math:`\bar{t}` 和无量纲 Green 函数 :math:`G_{ij}` ，
      表示 :math:`j` 方向的力激发的 :math:`i` 方向的位移。

      .. literalinclude:: run_lamb3/head_lamb3
         :language: text

      通过指定 **-S**，可输出更多内容：

      + 相对于源点坐标的一阶空间偏导 :math:`G_{ij,k'}`，即 :math:`\dfrac{\partial G_{ij}}{\partial x'_{k'}}`

      .. literalinclude:: run_lamb3/head_lamb3_source
         :language: text

      + 相对于接收点坐标的一阶空间偏导 :math:`G_{ij,k}`，即 :math:`\dfrac{\partial G_{ij}}{\partial x_k}`

      .. literalinclude:: run_lamb3/head_lamb3_receiver
         :language: text

      + 接收点和源点坐标的二阶混合偏导 :math:`G_{ij,k,k'}`，即
        :math:`\dfrac{\partial^2 G_{ij}}{\partial x_k \partial x'_{k'}}`

      .. literalinclude:: run_lamb3/head_lamb3_mixed
         :language: text

   .. group-tab:: Python

      Python 提供了函数 :func:`lamb3() <pygrt.utils.lamb3>`，直接返回四个数组
      ``G, Gs, Gr, Grs``，分别表示 Green 函数、源点一阶空间偏导、接收点一阶空间偏导
      和混合二阶空间偏导

      .. literalinclude:: run_lamb3/lamb3_plot_time.py
         :language: python
         :start-after: BEGIN LAMB3
         :end-before: END LAMB3

书中结果复现
--------------------

:download:`lamb3_plot_time.py <run_lamb3/lamb3_plot_time.py>`

.. figure:: run_lamb3/lamb3_0.1.svg
   :align: center

   复现了原书中的图 8.4.10

-----------

.. figure:: run_lamb3/lamb3_1.0.svg
   :align: center

   复现了原书中的图 8.4.11

-----------

.. figure:: run_lamb3/lamb3_5.0.svg
   :align: center

   复现了原书中的图 8.4.12

-----------

.. figure:: run_lamb3/lamb3_d1_0.1.svg
   :align: center

   复现了原书中的图 8.4.16

-----------

.. figure:: run_lamb3/lamb3_d2_0.1.svg
   :align: center

   复现了原书中的图 8.4.17

-----------

.. figure:: run_lamb3/lamb3_d3_0.1.svg
   :align: center

   复现了原书中的图 8.4.18

-----------

.. figure:: run_lamb3/lamb3_d1_1.0.svg
   :align: center

   复现了原书中的图 8.4.19

-----------

.. figure:: run_lamb3/lamb3_d2_1.0.svg
   :align: center

   复现了原书中的图 8.4.20

-----------

.. figure:: run_lamb3/lamb3_d3_1.0.svg
   :align: center

   复现了原书中的图 8.4.21

-----------

.. figure:: run_lamb3/lamb3_d1_5.0.svg
   :align: center

   复现了原书中的图 8.4.22

-----------

.. figure:: run_lamb3/lamb3_d2_5.0.svg
   :align: center

   复现了原书中的图 8.4.23

-----------

.. figure:: run_lamb3/lamb3_d3_5.0.svg
   :align: center

   复现了原书中的图 8.4.24

频域解和时域解的对比
-------------------------------

将第三类 Lamb 问题的格林函数对源点坐标的一阶偏导结果组合起来可以对应到频域解中的 EX, DD, DS, SS
震源结果，以下为对比验证（均已卷积阶跃函数）。

:download:`lamb3_plot_freq_time_source.py <run_lamb3/lamb3_plot_freq_time_source.py>`

.. figure:: run_lamb3/lamb3_compare_freq_time_source.svg
   :align: center

-----------

由于 PyGRT 的频域解中可以计算格林函数相对于接收点坐标的空间偏导，
因此这里我们可以进行更多的对比。在以下对比图中发现，对于频域解卷积了阶跃函数之后，
格林函数的 Gibbs 效应少了很多，然而对于格林函数的空间导数还是很明显，
这是因为经过理论推导，空间导数项转为了时间导数项。
如果想要更清晰的对比，可以自行对空间导数项再做一次积分然后绘制。

:download:`lamb3_plot_freq_time.py <run_lamb3/lamb3_plot_freq_time.py>`

.. figure:: run_lamb3/lamb3_compare_freq_time.svg
   :align: center

-----------

.. figure:: run_lamb3/lamb3_compare_freq_time_z.svg
   :align: center

-----------

.. figure:: run_lamb3/lamb3_compare_freq_time_r.svg
   :align: center

-----------

按照书中的思路，可进一步得到第三类 Lamb 问题的格林函数二阶混合偏导，
它可以组合为频域解中 EX、DD、DS、SS 震源对应的位移导数。下面给出对比验证；
为减弱频域结果的 Gibbs 振荡，对比脚本增加了一次时间积分，并去除了低频漂移。

:download:`lamb3_plot_freq_time_mixed.py <run_lamb3/lamb3_plot_freq_time_mixed.py>`

.. figure:: run_lamb3/lamb3_compare_freq_time_mixed_r.svg
   :align: center

-----------

.. figure:: run_lamb3/lamb3_compare_freq_time_mixed_z.svg
   :align: center

-----------
