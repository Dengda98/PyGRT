:author: 朱邓达
:date: 2026-08-20

应用在有限断层震源
======================

PyGRT 静态解的合成阶段目前已支持传入
`Coulomb 程序 <https://pubs.usgs.gov/of/2011/1060/>`_ 格式的有限断层。
此时震源不再是一个点源，PyGRT 会对断层面进行划分，
得到多个子断层，并将每个子断层视为在中心点的点源，
然后叠加各个子断层的合成结果，得到有限断层震源的合成结果。

有限断层文件
---------------
以下使用一个垂直走滑断层作为示例，文件 ``fault.inp`` 的内容如下：

.. literalinclude:: run_src_fault/fault.inp
    :language: text

前两行是注释，程序读取时会固定跳过。
关于格式的细节，详见 :doc:`/Module/static_syn` 模块和
:doc:`/Module/okada` 模块。

有限断层的空间展布大致如下
（:download:`plot_fault.py <run_src_fault/plot_fault.py>`）：

.. figure:: run_src_fault/fault.svg


计算静态格林函数库
-------------------------
显然，计算有限断层的结果至少需要 **多震源深度+多震中距** 的格林函数库。
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

计算有限断层激发的位移
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

:download:`Shell Scripts <run_src_fault/run.sh>` | :download:`Python Scripts <run_src_fault/run.py>`

.. literalinclude:: run_src_fault/run.py
    :language: python
    :start-after: BEGIN PLOT
    :end-before: END PLOT

.. figure:: run_src_fault/disp.svg
    :align: center

    颜色表示垂直位移 Z，箭头表示水平位移 E、N，粗黑线表示断层顶边走向

以上合成中使用 **-e** (C) 和 **calc_upar=True** (Python)，所以输出的 nc 文件中也包含位移偏导数，
而输出格式与之前的点源的情况没什么不同，因此应力等物理量的计算方式而之前的介绍完全一致，这里不再重复。
