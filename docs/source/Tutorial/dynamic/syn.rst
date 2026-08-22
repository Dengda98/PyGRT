:author: 朱邓达
:date: 2025-04-17

合成动态位移
=================

.. _warning_src_zdir:
.. warning:: 

    **震源机制参数中（如单力源、矩张量源）Z轴取向下为正。**

Python中合成动态位移的主函数为 :meth:`compute_syn() <pygrt.pymod.PyModel1D.compute_syn>` ，C模块为 :doc:`/Module/syn`。

使用上节计算的格林函数，合成动态位移（理论地震图）。方便起见，这里统一使用milrow模型，震源深度2km，场点位于地表，震中距10km的格林函数，方位角30°。
若仅使用已算好的格林函数，构造 :class:`~pygrt.pymod.PyModel1D` 时只需指定 ``grn``，无需再传 ``modelpath``：

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN REUSE GRN
    :end-before: END REUSE GRN

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

CLI 和 Python 函数中根据设置的不同震源参数自动推断震源类型。
一次只能设置一组震源专用参数，不完整或混用参数会报错。

脚本下载： :download:`Shell Scripts <run/run.sh>` | :download:`Python Scripts <run/run.py>`

.. tabs::

    .. group-tab:: C

        C中指定格林函数路径有以下两种方式：

        + 直接将 **-G** 指向震中距10km的格林函数子目录，例如
          ``-GGRN/milrow_2_0_10``。此时子目录已经确定了震源深度、台站深度和震中距，
          不能再设置 **-Ds/-Dr/-R**。
        + 将 **-G** 指向格林函数根目录，例如 ``-GGRN``，再使用必要的选项精确选择格林函数。
          本例中震源深度和台站深度各只有一个值，因此 **-Ds/-Dr** 可以省略；这里显式设置
          **-R10** 选择震中距为10km的子目录。

        本节的C示例使用第一种方式，脚本中也以注释形式给出了第二种方式。

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN SOURCES
            :end-before: END SYN SOURCES


    .. group-tab:: Python

        Python示例使用上节构造的格林函数根目录，并根据 ``dist`` 查找对应的格林函数；
        当根目录的对应维度只有一个值时，相应的选择参数可以省略，也可以显式设置正确的值。

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN SOURCES
            :end-before: END SYN SOURCES


.. figure:: run/syn_sources.svg
    :align: center
     
    不同震源参数对应的 Z 分量位移
     



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

        生成的时间函数会以SAC格式保存在对应路径中，文件名为 :file:`sig.sac`。
        其它时间函数以及具体参数用法详见 :meth:`compute_syn() <pygrt.pymod.PyModel1D.compute_syn>` 的 ``time_function`` 参数。

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

        Python 示例在合成后使用 :func:`stream_integral() <pygrt.utils.stream_integral>` /
        :func:`stream_diff() <pygrt.utils.stream_diff>` 做积分与微分。
        若希望在合成阶段完成，也可传入 ``integrate_order`` / ``differentiate_order``
        （分别对应 CLI 的 ``-I`` / ``-J``）。

.. figure:: run/syn_mt_intdif_Z.svg
   :align: center

.. figure:: run/syn_mt_intdif_R.svg
   :align: center

.. figure:: run/syn_mt_intdif_T.svg
   :align: center
