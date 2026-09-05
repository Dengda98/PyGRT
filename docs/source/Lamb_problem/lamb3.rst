:author: 朱邓达
:date: 2026-09-01

第三类 Lamb 问题
===================

第三类 Lamb 问题中，源点和观测点都位于半空间内部。

.. tabs::

   .. group-tab:: CLI

      C 程序 :command:`grt` 提供了模块 :doc:`/Module/lamb3` 求解第三类 Lamb 问题。

      .. literalinclude:: run_lamb3/run.sh
         :language: bash
         :start-after: BEGIN LAMB3
         :end-before: END LAMB3

      使用重定向将格林函数结果保存到文件 *lamb3.txt* 中，其内容格式类似于

      .. literalinclude:: run_lamb3/head_lamb3
         :language: text

      通过指定 **-S**，可将格林函数相对于源点/接收点坐标的一阶空间偏导分别写入对应路径。
      其内容格式类似于：

      + 相对于源点坐标

      .. literalinclude:: run_lamb3/head_lamb3_source
         :language: text

      + 相对于接收点坐标

      .. literalinclude:: run_lamb3/head_lamb3_receiver
         :language: text

   .. group-tab:: Python

      Python 提供了函数 :func:`lamb3() <pygrt.utils.lamb3>`，一次返回 Green 函数、
      相对于源点坐标的一阶空间偏导和相对于接收点坐标的一阶空间偏导。

      .. literalinclude:: run_lamb3/lamb3_plot_time.py
         :language: python
         :start-after: BEGIN LAMB3
         :end-before: END LAMB3

最后绘制计算得到的格林函数以及相对源点坐标的一阶空间导数。

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


关于格林函数对接收点坐标的一阶空间偏导的数学推导
---------------------------------------------------

以下沿用书中第 8 章的记号。上标 :math:`R` 和 :math:`S` 分别表示接收点坐标导数和
源点坐标导数。每个 :math:`F^R` 都与书中对应的源点导数积分同形，只需替换积分分子

接收点导数的积分组合
~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{aligned}
   \bar G^H_{ij,k}(\bar t)
   =\frac{\partial}{\partial\bar t}\Bigg\{\frac{1}{2}\Big[
   &\frac{\pi}{4}\left(F^R_{P,ij,k}+F^R_{S,ij,k}\right)\\
   &+\varsigma\left(F^R_{PP,ij,k}+F^R_{SS1,ij,k}
   -H(\theta'-\theta_c)\left(F^R_{SS2,ij,k}+F^R_{sPs,ij,k}\right)\right)
   -4\left(F^R_{PS,ij,k}+F^R_{SP,ij,k}\right)\Big]\Bigg\},
   \qquad k=1,2,3.
   \end{aligned}

直达 P、S 项
~~~~~~~~~~~~~

.. math::

   \boldsymbol D^{R,\lambda}_{,1}=-\boldsymbol D^{\lambda}_{,1'},
   \qquad
   \boldsymbol D^{R,\lambda}_{,2}=-\boldsymbol D^{\lambda}_{,2'},
   \qquad \lambda\in\{P,S\},

.. math::

   \boldsymbol D^{R,\lambda}_{,3}(x,x')
   =\left[\boldsymbol D^{\lambda}_{,3'}(x',x;\phi+\pi)\right]^T.

反射 PP、SS 项
~~~~~~~~~~~~~~~~

对 :math:`\lambda\in\{PP,SS\}`，源点竖直导数和接收点导数分别为

.. math::

   \boldsymbol M^{\lambda(\xi)}_{,3'}=-x\boldsymbol M^{\lambda(\xi)},
   \qquad
   \boldsymbol M^{R,\lambda(\xi)}_{,3}
   =\left[-x\boldsymbol M^{\lambda(\xi)}(\phi+\pi)\right]^T,
   \qquad \xi=1,2,

.. math::

   \boldsymbol M^{R,\lambda(\xi)}_{,1}=-\boldsymbol M^{\lambda(\xi)}_{,1'},
   \qquad
   \boldsymbol M^{R,\lambda(\xi)}_{,2}=-\boldsymbol M^{\lambda(\xi)}_{,2'}.

上式中的 :math:`\xi=1,2` 就是书中分母拆分后的两组分子；SS 这组分子同时用于
SS1、SS2 和 sPs

相应的两组积分核为

.. math::

   \begin{aligned}
   \mathrm{PP}:\quad
   &\frac{M^{R,PP(1)}_{ij,k}}{R'^{PP}\sqrt{Q_{PP}}}
   +\frac{M^{R,PP(2)}_{ij,k}}{R'^{PP}\sqrt{Q_{PP}}\eta_\beta},\\
   \mathrm{SS}:\quad
   &\frac{M^{R,SS(1)}_{ij,k}}{R'^{SS}\sqrt{Q_{SS}}}
   +\frac{M^{R,SS(2)}_{ij,k}}{R'^{SS}\sqrt{Q_{SS}}\eta_\alpha}.
   \end{aligned}

PS、SP 项
~~~~~~~~~~~

水平方向仍为

.. math::

   \boldsymbol M^{R,PS}_{,1}=-\boldsymbol M^{PS}_{,1'},
   \qquad
   \boldsymbol M^{R,PS}_{,2}=-\boldsymbol M^{PS}_{,2'},

.. math::

   \boldsymbol M^{R,SP}_{,1}=-\boldsymbol M^{SP}_{,1'},
   \qquad
   \boldsymbol M^{R,SP}_{,2}=-\boldsymbol M^{SP}_{,2'}.

竖直方向交换源点和接收点深度，并交换 PS 与 SP

.. math::

   \boldsymbol M^{R,PS}_{,3}(z_S,z_P;\phi)
   =\left[\boldsymbol M^{SP}_{,3'}(z_P,z_S;\phi+\pi)\right]^T,

.. math::

   \boldsymbol M^{R,SP}_{,3}(z_S,z_P;\phi)
   =\left[\boldsymbol M^{PS}_{,3'}(z_P,z_S;\phi+\pi)\right]^T.

源点竖直分子使用

.. math::

   \boldsymbol M^{PS}_{,3'}=-(x^2-1)\boldsymbol M^{PS},
   \qquad
   \boldsymbol M^{SP}_{,3'}=-(x^2+1)\boldsymbol M^{SP}.

因此 PS、SP 接收点竖直导数在书中转换项积分中直接替换为

.. math::

   F^{R,\lambda}_{ij,k}(\bar t)
   =\frac{k'^3r}{32R}\operatorname{Im}
   \int
   \frac{Q^{R,\lambda}_{ij,k}(x)}{x^4R_{PS}(x)\sqrt{W(x)}}\,\mathrm dx,
   \qquad \lambda\in\{PS,SP\},

.. math::

   Q^{R,\lambda}_{ij,k}(x)
   =(x^4-1)\mathcal F(x)M^{R,\lambda}_{ij,k}(x).
