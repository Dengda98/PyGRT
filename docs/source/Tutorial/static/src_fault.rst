:author: 朱邓达
:date: 2026-08-20

应用在有限断层震源
======================

PyGRT 静态解的合成阶段目前已支持传入
`Coulomb 程序 <https://pubs.usgs.gov/of/2011/1060/>`_ 格式的有限断层震源。
此时震源不再是一个点源，PyGRT 会对断层面进行划分，
得到多个子断层，并将每个子断层视为在中心点的点源，
然后叠加各个子断层的合成结果，得到有限断层震源的合成结果。

:download:`Shell Scripts <run_src_fault/run.sh>` | :download:`Python Scripts <run_src_fault/run.py>`

有限断层文件
---------------
以下使用一个垂直走滑断层作为示例，文件 ``fault.inp`` 的内容如下：

.. literalinclude:: run_src_fault/fault.inp
    :language: text

前两行是 Coulomb 表头：第一行首个 token 必须为 ``#``，其后按顺序给出 10 个字段标签，
第二行给出 11 个占位字段；常见表头中的 ``dip angle`` 可以是两个空白分隔的 token。
第 7 列表头精确写为 ``rake`` 时，第 7、8 列才按 rake/net slip 解释；文件名后缀不参与判断。
关于格式的细节，详见 :doc:`/Module/static_syn` 模块、
:doc:`/Module/okada` 模块以及 `Coulomb 程序 <https://pubs.usgs.gov/of/2011/1060/>`_ 程序使用手册。

有限断层的空间展布大致如下
（:download:`plot_fault.py <run_src_fault/plot_fault.py>`）：

.. figure:: run_src_fault/fault.svg
    :align: center


计算静态格林函数库
-------------------------
显然，计算有限断层震源的结果至少需要 **多震源深度+多震中距** 的格林函数库。
如果后续合成中不涉及多接收深度（假设接收深度都是 0.0），
则可以如下计算格林函数库。注意深度范围和震中距范围要覆盖你的研究区域。

.. tabs::

    .. group-tab:: C

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

    .. group-tab:: C

        .. literalinclude:: run_src_fault/run.sh
            :language: bash
            :start-after: BEGIN SYN
            :end-before: END SYN

        **注意， -X 指定北向坐标， -Y 指定东向坐标。**
        它们和 Coulomb 文件中 X/Y 的命名约定正好相反。

    .. group-tab:: Python

        .. literalinclude:: run_src_fault/run.py
            :language: python
            :start-after: BEGIN SYN
            :end-before: END SYN


.. literalinclude:: run_src_fault/run.py
    :language: python
    :start-after: BEGIN PLOT
    :end-before: END PLOT

.. figure:: run_src_fault/disp.svg
    :align: center

    颜色表示垂直位移 Z，箭头表示水平位移 E、N，粗黑线表示断层顶边走向

计算库伦应力
--------------------
以上合成中使用 **-e** (C) 和 **calc_upar=True** (Python)，所以输出的 nc 文件中也包含位移偏导数，
而输出格式与之前的点源的情况没什么不同，因此应力等物理量的计算方式而之前的介绍完全一致，这里不再重复。

得到计算得到位移偏导数后，就可以计算应力张量 -> 指定接收断层形态对应力张量进行投影 -> 计算库伦应力。

.. tabs::

    .. group-tab:: C

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
        :func:`compute_sproj() <pygrt.utils.compute_sproj>` 和
        :func:`compute_coulomb() <pygrt.utils.compute_coulomb>`。

.. figure:: run_src_fault/coulomb.svg
    :align: center
