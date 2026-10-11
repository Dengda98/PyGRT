:author: 朱邓达
:date: 2026-08-20

应用在有限断层震源
======================

静态合成通过 |sCoulomb| 格式的有限断层文件指定有限断层震源，
剖分后叠加各子断层中心点源的结果。

:download:`Shell Scripts <run_src_fault/run.sh>` | :download:`Python Scripts <run_src_fault/run.py>`

有限断层文件
---------------
以下使用一个垂直走滑断层作为示例，文件 ``fault.inp`` 的内容如下：

.. literalinclude:: run_src_fault/fault.inp
    :language: text

有限断层的空间展布大致如下
（:download:`plot_fault.py <run_src_fault/plot_fault.py>`）：

.. figure:: run_src_fault/fault.svg
    :align: center


计算静态格林函数库
-------------------------
一般使用 **多震源深度+多震中距** 的格林函数库。
如果后续合成中不涉及多接收深度（假设接收深度都是 0.0），
则可以如下计算格林函数库。注意深度范围和震中距范围要覆盖你的研究区域。

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run_src_fault/run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

    .. group-tab:: Python

        .. literalinclude:: run_src_fault/run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

计算位移及其空间偏导
------------------------------
在合成时， CLI 和 Python 函数中都有对应传入有限断层文件的参数，
此时点源相关的参数就不可再设置。
程序会根据格林函数库内各个参数的间隔最小值自动确定对有限断层的剖分，
也支持手动设置。

以下示例以二维接收网格为例，计算每个点上的位移以及位移偏导数。

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run_src_fault/run.sh
            :language: bash
            :start-after: BEGIN SYN
            :end-before: END SYN

        **注意， -X 指定北向坐标， -Y 指定东向坐标。**
        它们和断层文件中 X/Y 的命名约定正好相反，见 :ref:`finite-fault-file`。

    .. group-tab:: Python

        .. literalinclude:: run_src_fault/run.py
            :language: python
            :start-after: BEGIN SYN
            :end-before: END SYN

.. figure:: run_src_fault/disp.svg
    :align: center

    颜色表示垂直位移 Z，箭头表示水平位移 E、N，粗黑线表示断层顶边走向

计算库伦应力
--------------------
以上合成中使用 **-e** (CLI) 和 **calc_upar=True** (Python)，所以输出的 nc 文件中也包含位移偏导数，
而输出格式与之前的点源的情况没什么不同，因此应力等物理量的计算方式而之前的介绍完全一致，这里不再重复。

得到计算得到位移偏导数后，就可以计算应力张量 -> 指定接收断层形态对应力张量进行投影 -> 计算库伦应力。

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run_src_fault/run.sh
            :language: bash
            :start-after: BEGIN COULOMB
            :end-before: END COULOMB

        两个模块的详细说明请参见
        :doc:`static_sproj </Module/static_sproj>` 和
        :doc:`static_coulomb </Module/static_coulomb>`。

    .. group-tab:: Python

        .. literalinclude:: run_src_fault/run.py
            :language: python
            :start-after: BEGIN COULOMB
            :end-before: END COULOMB

        两个函数的详细说明请参见
        :func:`static_sproj() <pygrt.utils.static_sproj>` 和
        :func:`static_coulomb() <pygrt.utils.static_coulomb>`。

.. figure:: run_src_fault/coulomb.svg
    :align: center
