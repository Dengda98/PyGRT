:author: 朱邓达
:date: 2026-08-17

构建动态格林函数库与应用有限断层
========================================

.. admonition:: 本页文档示例下载
    :class: tip

    :download:`dynlib.tar.gz <dynlib.tar.gz>` |
    解压后进入示例目录，运行 ``bash run.sh`` 即可复现本页示例。

固定模型后，将不同震源深度、接收深度和震中距的格林函数保存成库，
即可反复用于不同震源机制、有限断层和接收点的合成。

构建动态格林函数库
------------------------

本节先考虑台站位于地表的情况，接收深度固定为 0 km，
因此只在 **震源深度** 和 **震中距** 两个维度上采样。

**格林函数库** 须满足以下条件：

* 目录下必须包含模型文件
* 只包含同一模型的结果，包含所选震源深度和震中距采样的全部组合
* 所有节点使用相同的采样点数、采样时间间隔等时间参数
* 包含目标震源类型需要的基本格林函数；若需计算应力等物理量，还须包含空间导数

采样范围须覆盖所有子源的深度，以及各子源到地表接收点的实际震中距。
采样空间间隔将影响插值精度，但采样过多会要求更多的计算时间和磁盘存储空间，因此要适当选择。
各轴可以非等距。

以下示例计算震源深度 2～6 km、震中距 2～36 km 的所有组合，
时间采样为 512 点、间隔 0.1 s，同时计算格林函数及其空间导数。

.. tabs::

    .. group-tab:: CLI

        使用 :doc:`/Module/greenfn` 一次建库，**-Dr0** 固定地表接收点，
        **-Ds/-R** 指定震源深度和震中距采样，
        **-e** 同时计算空间导数。

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

    .. group-tab:: Python

        调用 :meth:`greenfn() <pygrt.pymod.PyModel1D.greenfn>`
        一次计算全部深度和距离组合。

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

输出根目录为 ``GRN``，每个深度与距离组合对应一个 SAC 子目录：

.. code-block:: text

    GRN/
    ├── mod3_2_0_2/
    │   ├── DDZ.sac
    │   └── ...
    ├── mod3_2_0_4/
    └── ...

子目录命名为 ``{model}_{depsrc}_{deprcv}_{dist}``，
分别对应模型名、震源深度、接收深度和震中距。
从库根目录合成时，默认线性插值，也支持最近邻查询。
目标深度和震中距不必与采样节点重合，但不能超出库的范围。

应用有限断层震源
------------------------------

**使用有限断层震源时，输出只能采用 ZNE 坐标系。**

有限断层的剖分、各子断层的时间函数和破裂延迟可以自行定义，也可以从其他来源获得后调整为 |sCoulomb| 格式。
本例用 :doc:`/Module/faultdiv` 模块辅助快速生成 1×1 km 的子断层及简单的破裂延迟，
各子断层采用上升、下降各 1 s 的余弦时间函数，破裂从断层面中心以 2.5 km/s 向外传播。
文件前几行如下：

.. literalinclude:: fault_rupture.inp
    :language: text
    :lines: 1-5


.. figure:: fault.svg
    :align: center

    剖分后震源断层的水平投影和三维展布

本例在地表台站计算阶跃位移及空间导数，
随后计算应力张量，投影到走向 90°、倾角 60°、滑动角 0° 的接收面，
得到法向应力和沿滑动方向的切向（剪）应力，最后以等效摩擦系数 0.4 计算库伦应力。

.. tabs::

    .. group-tab:: CLI

        使用 :doc:`/Module/syn` 的 **-C** 读取有限震源，**+i0/0** 保留文件中的现有剖分。
        不设置全局 **-D**，以保留逐行时间函数和延迟。
        库查询默认线性插值，**-i0** 改为最近邻。

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN SRC SYN
            :end-before: END SRC SYN

    .. group-tab:: Python

        通过 :meth:`syn() <pygrt.pymod.PyModel1D.syn>` 的 *src_fault* 读取有限震源，
        *src_fault_size=(0, 0)* 保留现有剖分。
        不设置全局 *time_function*，以保留逐行时间函数和延迟。
        库查询默认线性插值，*interpolate=False* 改为最近邻。

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN SRC SYN
            :end-before: END SRC SYN

.. figure:: syn_ff.svg
    :align: center

    三分量位移

.. figure:: coulomb.svg
    :align: center

    切向应力、法向应力和库伦应力变化

应用有限断层接收
------------------------------

:doc:`/Module/syn` 和 :meth:`syn() <pygrt.pymod.PyModel1D.syn>`
也支持传入有限接收断层， **本质上是批量定义多个接收点**。
可以在整条断层中心接收，也可剖分后在各子断层中心接收。

**上述格林函数库的接收深度固定为 0 km。由于接收断层包含地下接收点，
就必须在“接收点深度”维度上扩增格林函数库。**
同时仍需确保各子源到接收点的实际震中距位于库的范围内。

接收断层与震源断层使用相同的 |sCoulomb| 文件格式，但仅支持 **Kode=100**。
滑动量不参与位移计算，行末的震源时间函数和延迟会被忽略。
几何和滑动角随结果保存，供后续应力投影使用；
若滑动角未定义，须在投影阶段补充，也可用 **rake** 格式预先给定。

同一次合成也可同时使用有限震源和有限接收断层，当然这通常就要求格林函数库在三个空间维度上展布。
后处理方式与点源相同。
