:author: 朱邓达
:date: 2025-11-14

第一类 Lamb 问题
===================

第一类 Lamb 问题是指，在半空间模型中，源点和场点均位于地表，求解场点记录到的位移。

.. tabs::  

   .. group-tab:: C
      
      C 程序 :command:`grt` 提供了模块 :doc:`/Module/lamb1` 求解第一类 Lamb 问题。

      .. literalinclude:: run/run.sh
         :language: bash
         :start-after: BEGIN LAMB1
         :end-before: END LAMB1

      使用重定向将结果保存到文件 *lamb1.txt* 中，其内容格式类似于

      .. literalinclude:: run/head_lamb1
         :language: text

   .. group-tab:: Python    

      Python 提供了函数 :func:`solve_lamb1() <pygrt.utils.solve_lamb1>` 求解第一类 Lamb 问题。

      .. literalinclude:: run/lamb1_plot_time.py
         :language: python
         :start-after: BEGIN LAMB1
         :end-before: END LAMB1

最后绘制计算得到的格林函数。

:download:`lamb1_plot_time.py <run/lamb1_plot_time.py>`

.. figure:: run/lamb1_time.svg
   :align: center


频域解和时域解的对比
-------------------------------

对比观察可发现频域解在波形突变出有明显的 Gibbs 效应。

:download:`lamb1_plot_freq_time.py <run/lamb1_plot_freq_time.py>`

.. figure:: run/lamb1_compare_freq_time.svg
   :align: center


