:author: 朱邓达
:date: 2026-09-16

lamb 模块
===================

:doc:`/Module/lamb` 模块用于计算均匀各向同性弹性半空间中的动态全波解的解析解。
其底层由 :doc:`/Module/lamb1`、:doc:`/Module/lamb2` 和
:doc:`/Module/lamb3` 支持，通过分量整理得到各类震源激发的三分量位移及其空间偏导数的 SAC 记录，
输出的文件形式对标 :doc:`/Module/syn` 模块。

下面在 :math:`V_P=8.0` km/s、:math:`V_S=4.62` km/s、
:math:`\rho=3.3` g/cm\ :sup:`3` 的均匀半空间中，计算震源深度为 5 km 的剪切源 (*strike/dip/rake=33/50/120*)，
在深度为 1 km，水平距离为 15 km ，方位角为 30 ° 的接收点处记录的位移，
解析解和数值解都卷积宽度为 0.2 s 的三角波，并在绘图时积分一次得到阶跃型位移。

.. tabs::

   .. group-tab:: CLI

      具体用法可见 :doc:`/Module/lamb` 模块。

      .. literalinclude:: run_lamb/run.sh
         :language: bash
         :start-after: BEGIN LAMB
         :end-before: END LAMB

   .. group-tab:: Python

      具体用法可见 :func:`lamb() <pygrt.utils.lamb>` 函数。

      .. literalinclude:: run_lamb/plot_compare.py
         :language: python
         :start-after: BEGIN LAMB
         :end-before: END LAMB

:download:`plot_compare.py <run_lamb/plot_compare.py>`

.. figure:: run_lamb/lamb_compare.svg
   :align: center
   :width: 90%

   均匀半空间中 Lamb 解析解与 greenfn + syn 结果的 ZNE 分量对比
