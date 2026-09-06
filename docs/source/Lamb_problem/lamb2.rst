:author: 朱邓达
:date: 2026-08-31

第二类 Lamb 问题
===================

第二类 Lamb 问题是指，在半空间模型中，震源位于地下、接收点位于地表的情形。
不过地表源、地下接收可由互易定理得到。

.. tabs::

   .. group-tab:: CLI

      C 程序 :command:`grt` 提供了模块 :doc:`/Module/lamb2` 求解第二类 Lamb 问题。

      .. literalinclude:: run_lamb2/run.sh
         :language: bash
         :start-after: BEGIN LAMB2
         :end-before: END LAMB2

      使用重定向将格林函数结果保存到文件 *lamb2.txt* 中，其内容格式类似于

      .. literalinclude:: run_lamb2/head_lamb2
         :language: text

      通过指定 **-S**，可将格林函数相对于源点/接收点坐标的一阶空间偏导分别写入对应路径。
      其内容格式类似于：

      + 相对于源点坐标

      .. literalinclude:: run_lamb2/head_lamb2_source
         :language: text

      + 相对于接收点坐标

      .. literalinclude:: run_lamb2/head_lamb2_receiver
         :language: text


   .. group-tab:: Python

      Python 提供了函数 :func:`lamb2() <pygrt.utils.lamb2>`，一次返回 Green 函数、
      相对于源点坐标的一阶空间偏导和相对于接收点坐标的一阶空间偏导。

      .. literalinclude:: run_lamb2/lamb2_plot_time.py
         :language: python
         :start-after: BEGIN LAMB2
         :end-before: END LAMB2


最后绘制计算得到的格林函数以及相对源点坐标的一阶空间导数。

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


关于格林函数对接收点坐标的一阶空间偏导的数学推导
---------------------------------------------------

以下沿用书中第 7 章的记号。上标 :math:`R` 和 :math:`S` 分别表示接收点坐标导数和
源点坐标导数；P 项使用 :math:`M`，S1、S2 和 S-P 项使用 :math:`N`。每个
:math:`F^R` 都与书中对应的源点导数积分同形，只需替换积分分子

接收点导数的积分组合
~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{aligned}
   \bar G^H_{ij,k}(\bar t)
   =\frac{\partial}{\partial\bar t}\Big[&F^R_{P,ij,k}(\bar t)+F^R_{S1,ij,k}(\bar t)
    -H(\theta-\theta_c)\{F^R_{S2,ij,k}(\bar t)+F^R_{S-P,ij,k}(\bar t)\}\Big],
   \qquad k=1,2,3.
   \end{aligned}

矩阵关系
~~~~~~~~

对 :math:`\xi=1,2`，水平接收点导数直接由源点导数得到

.. math::

   M^{R,P(\xi)}_{,1}=-M^{P(\xi)}_{,1'},
   \qquad
   M^{R,P(\xi)}_{,2}=-M^{P(\xi)}_{,2'},

.. math::

   N^{R,S(\xi)}_{,1}=-N^{S(\xi)}_{,1'},
   \qquad
   N^{R,S(\xi)}_{,2}=-N^{S(\xi)}_{,2'}.

对 :math:`A=M^P` 或 :math:`N^S`，接收点竖直导数矩阵为

.. math::

   A^{R,(\xi)}_{,3}=
   \begin{bmatrix}
   A^{(\xi)}_{31,1'}&A^{(\xi)}_{32,1'}&A^{(\xi)}_{33,1'}\\
   A^{(\xi)}_{31,2'}&A^{(\xi)}_{32,2'}&A^{(\xi)}_{33,2'}\\
   A^{R,(\xi)}_{31,3}&A^{R,(\xi)}_{32,3}&A^{R,(\xi)}_{33,3}
   \end{bmatrix}.

未按书中 :math:`R'` 与 :math:`R'\eta` 拆分时，新增的最后一行为

.. math::

   \begin{aligned}
   \left(M^R_{31,3},M^R_{32,3},M^R_{33,3}\right)
   &=-2\chi(q^2-p^2)\left(q\eta_\beta c_\phi,\ q\eta_\beta s_\phi,\ \eta_\alpha\eta_\beta\right),\\
   \left(N^R_{31,3},N^R_{32,3},N^R_{33,3}\right)
   &=\chi\gamma_S(q^2-p^2)\left(-q\eta_\beta c_\phi,\ -q\eta_\beta s_\phi,\ 1\right).
   \end{aligned}

拆分后的积分核和新增分子
~~~~~~~~~~~~~~~~~~~~~~~~~~

P 项的两组积分分子出现在

.. math::

   \frac{M^{R,P(1)}_{ij,k}}{R'^P\sqrt{Q_P}}
   +\frac{M^{R,P(2)}_{ij,k}}{R'^P\sqrt{Q_P}\eta_\beta},

其中接收点竖直导数最后一行的新增分子为

.. math::

   \begin{aligned}
   \left(M^{R,P(1)}_{31,3},M^{R,P(1)}_{32,3},M^{R,P(1)}_{33,3}\right)
   &=-8\chi x^2g\eta_\beta^2(q^2-p^2)\left(qc_\phi,\ qs_\phi,\ x\right),\\
   \left(M^{R,P(2)}_{31,3},M^{R,P(2)}_{32,3},M^{R,P(2)}_{33,3}\right)
   &=-2\chi x\eta_\beta^2\gamma_P^2(q^2-p^2)\left(qc_\phi,\ qs_\phi,\ x\right).
   \end{aligned}

S 项以及 S-P 项的两组积分分子出现在

.. math::

   \frac{N^{R,S(1)}_{ij,k}}{R'^S\sqrt{Q_S}}
   +\frac{N^{R,S(2)}_{ij,k}}{R'^S\sqrt{Q_S}\eta_\alpha},

其中接收点竖直导数最后一行的新增分子为

.. math::

   \begin{aligned}
   \left(N^{R,S(1)}_{31,3},N^{R,S(1)}_{32,3},N^{R,S(1)}_{33,3}\right)
   &=\chi x\gamma_S^3\left(-xq c_\phi,\ -xq s_\phi,\ q^2-p^2\right),\\
   \left(N^{R,S(2)}_{31,3},N^{R,S(2)}_{32,3},N^{R,S(2)}_{33,3}\right)
   &=4\chi x^2g'\eta_\alpha^2\gamma_S\left(-xq c_\phi,\ -xq s_\phi,\ q^2-p^2\right).
   \end{aligned}

S1、S2 和 S-P 共用上述 :math:`N` 分子，区别仍只在书中规定的积分路径、积分上下限
和基本积分
