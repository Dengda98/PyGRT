:author: 朱邓达
:date: 2026-08-20

应用在有限断层接收
======================

与上一节 :doc:`src_fault` 类似，PyGRT 静态解的合成阶段目前还支持通过传入
`Coulomb 程序 <https://pubs.usgs.gov/of/2011/1060/>`_ 格式的有限断层来指定接收点。
类似的，PyGRT 会对有限断层进行划分，将每个子断层中心点位置视为接收点，计算位移及偏导数。

与作为震源的有限断层不一样，有限接收断层中和滑动相关的量在计算位移时并不会用到，
但在后续计算流程中（如库伦应力）会用到，因此合成阶段的输出文件会记录有限接收断层的形态。

以下示例以一个点源（剪切源）为例。为简单起见，假设震源深度 5 km，计算在一个垂直走滑断层上的位移分布
（这相当于做了一个纵剖面）。

:download:`Shell Scripts <run_rcv_fault/run.sh>` | :download:`Python Scripts <run_rcv_fault/run.py>`


有限断层文件
------------------
该有限断层文件内部为

.. literalinclude:: run_rcv_fault/rcv_fault.inp
    :language: text

计算静态格林函数库
-----------------------

由于我们当前仅考虑一个固定深度的点源，因此建立格林函数库要求 **多接收深度+多震中距** 。

.. tabs::

    .. group-tab:: C

        .. literalinclude:: run_rcv_fault/run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

    .. group-tab:: Python

        .. literalinclude:: run_rcv_fault/run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

计算有限断层上接收的位移
------------------------------
程序会根据格林函数库内各个参数的间隔最小值自动确定对有限断层的剖分，也支持手动设置。
输出文件中会记录每个子断层中心点的坐标以及对应的位移及其偏导数结果。

.. tabs::

    .. group-tab:: C

        .. literalinclude:: run_rcv_fault/run.sh
            :language: bash
            :start-after: BEGIN SYN
            :end-before: END SYN

    .. group-tab:: Python

        .. literalinclude:: run_rcv_fault/run.py
            :language: python
            :start-after: BEGIN SYN
            :end-before: END SYN

这个可视化对三维空间内的位移并不直观，这里只是一个示例（:download:`plot_rcv_fault.py <run_rcv_fault/plot_rcv_fault.py>`），
建议用更专业的工具来进行绘制。当然，实际研究中也通常不使用有限接收断层来看位移，而是用位移偏导数。

.. figure:: run_rcv_fault/rcv_fault_zne.svg
    :align: center

    在作为接收平面的断层面上的三分量位移


以上教程只是抛砖引玉。如果对于震源和接收点都是有限断层的情况，那么建立格林函数库就需要 **多震源深度+多接收深度+多震中距** ，
这也是最通用的格林函数库建立方式。感兴趣的研究者可以自行测试，也欢迎提供样例。
