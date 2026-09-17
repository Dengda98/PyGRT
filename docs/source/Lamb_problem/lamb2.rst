:author: 朱邓达
:date: 2026-08-31

第二类 Lamb 问题（lamb2 模块）
==================================

第二类 Lamb 问题是指，在半空间模型中，震源位于地下、接收点位于地表的情形。
不过地表源、地下接收可由互易定理得到。

.. tabs::

   .. group-tab:: CLI

      :command:`grt` 命令提供了模块 :doc:`/Module/lamb2` 求解第二类 Lamb 问题。

      .. literalinclude:: run_lamb2/run.sh
         :language: bash
         :start-after: BEGIN LAMB2
         :end-before: END LAMB2

      使用重定向将结果保存到文件 *lamb2.txt* 中，其内容格式如下，
      记录了无量纲时间 :math:`\bar{t}` 和无量纲 Green 函数 :math:`G_{ij}` ，
      表示 :math:`j` 方向的力激发的 :math:`i` 方向的位移。

      .. literalinclude:: run_lamb2/head_lamb2
         :language: text

      通过指定 **-S**，可输出更多内容：

      + 相对于源点坐标的一阶空间偏导 :math:`G_{ij,k'}`，即 :math:`\dfrac{\partial G_{ij}}{\partial x'_{k'}}`

      .. literalinclude:: run_lamb2/head_lamb2_source
         :language: text

      + 相对于接收点坐标的一阶空间偏导 :math:`G_{ij,k}`，即 :math:`\dfrac{\partial G_{ij}}{\partial x_k}`

      .. literalinclude:: run_lamb2/head_lamb2_receiver
         :language: text

      + 接收点和源点坐标的二阶混合偏导 :math:`G_{ij,k,k'}`，即
        :math:`\dfrac{\partial^2 G_{ij}}{\partial x_k \partial x'_{k'}}`

      .. literalinclude:: run_lamb2/head_lamb2_mixed
         :language: text

   .. group-tab:: Python

      Python 提供了函数 :func:`lamb2() <pygrt.utils.lamb2>`，直接返回四个数组
      ``G, Gs, Gr, Grs``，分别表示 Green 函数、源点一阶空间偏导、接收点一阶空间偏导
      和混合二阶空间偏导

      .. literalinclude:: run_lamb2/lamb2_plot_time.py
         :language: python
         :start-after: BEGIN LAMB2
         :end-before: END LAMB2

书中结果复现
--------------------

:download:`lamb2_plot_time.py <run_lamb2/lamb2_plot_time.py>`

.. figure:: run_lamb2/lamb2.svg
   :align: center

   复现了原书中的图 7.4.6

-----------

.. figure:: run_lamb2/lamb2_d1.svg
   :align: center

   复现了原书中的图 7.4.10

-----------

.. figure:: run_lamb2/lamb2_d2.svg
   :align: center

   复现了原书中的图 7.4.11

-----------

.. figure:: run_lamb2/lamb2_d3.svg
   :align: center

   复现了原书中的图 7.4.12

频域解和时域解的对比
-------------------------------

将第二类 Lamb 问题的格林函数对源点坐标的一阶偏导结果组合起来可以对应到频域解中的 EX, DD, DS, SS
震源结果，以下为对比验证（均已卷积阶跃函数）。

:download:`lamb2_plot_freq_time_source.py <run_lamb2/lamb2_plot_freq_time_source.py>`

.. figure:: run_lamb2/lamb2_compare_freq_time_source.svg
   :align: center

-----------

由于 PyGRT 的频域解中可以计算格林函数相对于接收点坐标的空间偏导，
因此这里我们可以进行更多的对比。在以下对比图中发现，对于频域解卷积了阶跃函数之后，
格林函数的 Gibbs 效应少了很多，然而对于格林函数的空间导数还是很明显，
这是因为经过理论推导，空间导数项转为了时间导数项。
如果想要更清晰的对比，可以自行对空间导数项再做一次积分然后绘制。

:download:`lamb2_plot_freq_time.py <run_lamb2/lamb2_plot_freq_time.py>`

.. figure:: run_lamb2/lamb2_compare_freq_time.svg
   :align: center

-----------

.. figure:: run_lamb2/lamb2_compare_freq_time_z.svg
   :align: center

-----------

.. figure:: run_lamb2/lamb2_compare_freq_time_r.svg
   :align: center

---------------

按照书中的思路，可进一步得到第二类 Lamb 问题的格林函数二阶混合偏导，
它可以组合为频域解中 EX、DD、DS、SS 震源对应的位移导数。下面给出对比验证；
为减弱频域结果的 Gibbs 振荡，对比脚本增加了一次时间积分，并去除了低频漂移。

:download:`lamb2_plot_freq_time_mixed.py <run_lamb2/lamb2_plot_freq_time_mixed.py>`

.. figure:: run_lamb2/lamb2_compare_freq_time_mixed_r.svg
   :align: center

-----------

.. figure:: run_lamb2/lamb2_compare_freq_time_mixed_z.svg
   :align: center

-----------
