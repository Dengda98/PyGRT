:author: 朱邓达
:date: 2025-04-17

合成动态位移
=================

.. admonition:: 本页文档示例下载
    :class: tip

    :download:`syn.tar.gz <syn.tar.gz>` |
    解压后进入示例目录，运行 ``bash run.sh`` 即可复现本页示例。

.. _warning_src_zdir:
.. warning:: 

    **震源机制参数中（如单力源、矩张量源）Z轴取向下为正。**

Python 中合成动态位移的主函数为 :meth:`syn() <pygrt.pymod.PyModel1D.syn>`，
命令行中使用 :doc:`/Module/syn` 模块。

使用上节计算的格林函数，合成动态位移（理论地震图）。方便起见，这里统一使用milrow模型，震源深度2km，场点位于地表，震中距10km的格林函数，方位角30°。
若仅使用已算好的格林函数，构造 :class:`~pygrt.pymod.PyModel1D` 时只需指定 ``grn``，无需再传 ``modelpath``：

.. literalinclude:: run.py
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

其中 :math:`D(t)` 为完整的震源时间函数，描述了力加载过程或断层破裂过程，:math:`*` 表示卷积，
:math:`A_m` 为与方位角和震源机制相关的方向因子，其中 :math:`u_z, u_r` 的方向因子相同，
而 :math:`u_\theta` 的方向因子满足

.. math:: 

    A_{m+3} = \frac{d A_m}{d (m\theta)}, m=1,2
    
其中 :math:`m` 为阶数，:math:`\theta` 为方位角。


.. note:: 

    合成位移的结果单位为 :math:`\text{cm}`。


.. _dynamic-syn-sources:

不同震源
--------------

CLI 和 Python 函数中根据设置的不同震源参数自动推断震源类型。
一次只能设置一组震源专用参数，不完整或混用参数会报错。
本节示例均为点源。爆炸源、剪切源、张裂源和矩张量源采用阶跃响应，单力源采用脉冲响应。


.. tabs::

    .. group-tab:: CLI

        :command:`grt` 命令在命令行指定格林函数路径有以下两种方式：

        + 直接将 **-G** 指向震中距10km的格林函数子目录，例如
          ``-GGRN/milrow_2_0_10``。此时子目录已经确定了震源深度、台站深度和震中距，
          不能再设置 **-Ds/-Dr/-R**。
        + 将 **-G** 指向库根目录，例如 ``-GGRN``，给出目标深度和震中距，
          默认使用线性插值，也可用 **-i0** 改为最近邻查询。
          线性插值先将角点波形的 P 波初至按整数采样点对齐到目标位置的初至，再加权叠加。
          目标值不必与采样节点重合，但须在库的范围内。
          本例震源深度和台站深度各只有一个值，**-Ds/-Dr** 可以省略；
          例如 **-R7** 可在已有的 5 km 和 8 km 节点之间合成 7 km 处的位移。

        本节的C示例使用第一种方式，脚本中也以注释形式给出了第二种方式。

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN SYN SOURCES
            :end-before: END SYN SOURCES


    .. group-tab:: Python

        Python 示例从库根目录按目标深度和震中距查询，默认线性插值，
        设置 ``interpolate=False`` 可改为最近邻。
        例如 ``dist=7.0`` 不需要已有 7 km 节点，只需位于库的距离范围内。
        只有一个采样值的维度可以省略对应参数，该维度仍只适用于这个值。

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN SYN SOURCES
            :end-before: END SYN SOURCES


.. figure:: syn_sources.svg
    :align: center
     
    不同震源参数对应的 Z 分量位移
     


分量旋转
---------------------
**PyGRT** 计算默认输出为ZRT分量（柱坐标系），可以设置参数以输出ZNE分量，这里以剪切源为例，

.. tabs::  

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN ZNE
            :end-before: END ZNE

    .. group-tab:: Python 

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN ZNE
            :end-before: END ZNE


.. figure:: syn_dc_zne.svg
   :align: center


卷积时间函数
---------------------
**PyGRT** 内置抛物波、梯形波、非对称余弦波等震源时间函数，也支持自定义时间函数。
这里以剪切源为例，计算阶跃响应，卷积面积归一化的时间函数，
其累积描述断层错动过程，合成结果为位移。
自定义时间函数支持一列振幅，或两列时间（s）、振幅。
两列输入会线性插值到采样间隔 *dt*。程序会对自定义时间函数自动进行面积归一化。

下面对上述 **-M33/50/120** 剪切源分别卷积三角波和非对称余弦波，两者的持续时间均为 0.6 s。
三角波的上升段和下降段各为 0.3 s，平台段为零；非对称余弦波的上升段为 0.1 s，
下降段为 0.5 s，峰值出现在 0.1 s。

.. tabs::  

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN TIME FUNC
            :end-before: END TIME FUNC

        生成的时间函数会以SAC格式保存在对应路径中，文件名为 :file:`sig.sac`。

    .. group-tab:: Python 

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN TIME FUNC
            :end-before: END TIME FUNC

        生成的时间函数会以SAC格式保存在对应路径中，文件名为 :file:`sig.sac`。

.. figure:: syn_dc_trig.svg
   :align: center

   三角波矩率形状及其合成位移（上升 0.3 s，平台 0 s，下降 0.3 s）

.. figure:: syn_dc_cosine.svg
   :align: center

   非对称余弦波矩率形状及其合成位移（上升 0.1 s，下降 0.5 s）


位移对时间积分、微分
--------------------------------
合成时可以指定时间积分、微分次数，这些操作已经在 :ref:`不同震源 <dynamic-syn-sources>` 中使用。
如果疑惑积分/微分操作对结果物理单位的改变，详见 :doc:`../disp_vel/disp_vel`。


.. tabs::  

    .. group-tab:: CLI

        :doc:`/Module/syn` 模块的 **-I** 和 **-J** 分别指定时间积分和微分次数，
        例如 **-I1** 表示积分一次。

    .. group-tab:: Python 

        :meth:`syn() <pygrt.pymod.PyModel1D.syn>` 的 ``integrate_order`` 和 ``differentiate_order``
        分别指定时间积分和微分次数。
