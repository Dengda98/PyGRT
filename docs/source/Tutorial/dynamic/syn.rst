:author: 朱邓达
:date: 2025-04-17

合成动态位移
=================

.. _warning_src_zdir:
.. warning:: 

    **震源机制参数中（如单力源、矩张量源）Z轴取向下为正。**

Python中合成动态位移的主函数为 :func:`gen_syn_from_gf_*() <pygrt.utils.gen_syn_from_gf_DC>` （\*表示对不同震源）  ，C模块为 :doc:`/Module/syn`。

使用上节计算的格林函数，合成动态位移（理论地震图）。方便起见，这里统一使用milrow模型，震源深度2km，场点位于地表，震中距10km的格林函数，方位角30°。

在已知三分量格林函数 :math:`W_m(t), Q_m(t), V_m(t)` 后，合成三分量位移 :math:`u_z(t), u_r(t), u_\theta (t)` 的公式为

.. math:: 

   \left\{
   \begin{aligned}
    u_z(t) &= D(t) * \left[ \sum_{m=0}^{m=2} A_m W_m(t) \right] \\
    u_r(t) &= D(t) * \left[ \sum_{m=0}^{m=2} A_m Q_m(t) \right] \\
    u_\theta (t) &= D(t) * \left[ \sum_{m=1}^{m=2} A_{m+3} V_m(t) \right]
    \end{aligned}
   \right.

其中 :math:`D(t)` 为震源时间函数，:math:`*` 表示卷积，:math:`A_m` 为与方位角和震源机制相关的方向因子，其中 :math:`u_z, u_r` 的方向因子相同，而 :math:`u_\theta` 的方向因子满足 

.. math:: 

    A_{m+3} = \frac{d A_m}{d (m\theta)}, m=1,2
    
其中 :math:`m` 为阶数，:math:`\theta` 为方位角。


.. note:: 

    合成位移的结果单位为 :math:`\text{cm}`。


不同震源
--------------

以下绘图使用Python绘制，绘图函数如下：

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN plot func
    :end-before: END plot func


爆炸源
~~~~~~~~~~~~~~~~~
标量矩 1e24 dyne·cm。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN EX
            :end-before: END SYN EX

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN EX
            :end-before: END SYN EX


.. figure:: run/syn_ex.svg
   :align: center



单力源
~~~~~~~~~~~~~~~~~
北向力 :math:`f_N=1`，东向力 :math:`f_E=-0.5`，垂直向下的力 :math:`f_Z=2`，单位 1e16 dyne。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN SF
            :end-before: END SYN SF

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN SF
            :end-before: END SYN SF


.. figure:: run/syn_sf.svg
   :align: center


剪切源
~~~~~~~~~~~~~~
断层走向33°，倾角50°，滑动角120°，标量矩 1e24 dyne·cm。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN DC
            :end-before: END SYN DC

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN DC
            :end-before: END SYN DC


.. figure:: run/syn_dc.svg
   :align: center


张位错源
~~~~~~~~~~~~~~
断层走向33°，倾角50°，标量矩 1e24 dyne·cm。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN TS
            :end-before: END SYN TS

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN TS
            :end-before: END SYN TS


.. figure:: run/syn_ts.svg
   :align: center


矩张量源
~~~~~~~~~~~~~~
:math:`M_{xx}=0.1, M_{xy}=-0.2, M_{xz}=1.0, M_{yy}=0.3, M_{yz}=-0.5, M_{zz}=-2.0`，单位 1e24 dyne·cm， **其中X为北向，Y为东向，Z为垂直向下**。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN MT
            :end-before: END SYN MT

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN MT
            :end-before: END SYN MT


.. figure:: run/syn_mt.svg
   :align: center



分量旋转
---------------------
**PyGRT** 计算默认输出为ZRT分量（柱坐标系），可以设置参数以输出ZNE分量，这里以剪切源为例，

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN ZNE
            :end-before: END ZNE

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN ZNE
            :end-before: END ZNE


.. figure:: run/syn_dc_zne.svg
   :align: center



卷积时间函数
---------------------
**PyGRT** 内置了一些震源时间函数，例如抛物波、梯形波、雷克子波或自定义，这里以单力源为例。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN TIME FUNC
            :end-before: END TIME FUNC

        生成的时间函数会以SAC格式保存在对应路径中，文件名为 :file:`sig.sac`。 其它时间函数以及具体参数用法详见 :doc:`/Module/syn` 模块。

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN TIME FUNC
            :end-before: END TIME FUNC

        其它时间函数以及具体参数用法可在 :py:mod:`pygrt.signals` 模块中查看函数参数。

.. figure:: run/syn_sf_trig.svg
   :align: center



位移对时间积分、微分
--------------------------------
这里以矩张量源为例。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN INT DIF
            :end-before: END INT DIF


    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN INT DIF
            :end-before: END INT DIF


.. figure:: run/syn_mt_intdif_Z.svg
   :align: center

.. figure:: run/syn_mt_intdif_R.svg
   :align: center

.. figure:: run/syn_mt_intdif_T.svg
   :align: center