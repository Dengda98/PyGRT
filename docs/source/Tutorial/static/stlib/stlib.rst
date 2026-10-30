:author: 朱邓达
:date: 2026-08-17

构建静态格林函数库与应用有限断层
========================================

.. admonition:: 本页文档示例下载
    :class: tip

    :download:`stlib.tar.gz <stlib.tar.gz>` |
    解压后进入示例目录，运行 ``bash run.sh`` 即可复现本页示例。

静态格林函数库同样可以反复用于不同震源机制和接收位置的合成，
所有采样结果保存在一个 |NetCDF| 文件中。

构建静态格林函数库
------------------------

本节先考虑台站位于地表的情况，接收深度固定为 0 km，
因此只在 **震源深度** 和 **震中距** 两个维度上采样。

使用同一模型，计算所选震源深度和震中距采样的全部组合。
库的范围须覆盖所有子源的深度，以及各子源到地表接收点的实际震中距；
采样间隔应足够小，以保证插值精度。

本例计算震源深度 2～6 km、震中距 0～36 km 的所有组合，
同时计算位移空间导数，供后续计算应力等物理量。

.. tabs::

    .. group-tab:: CLI

        使用 :doc:`/Module/static_greenfn` 模块一次建库，**-Dr0** 固定地表接收点，
        **-Ds/-R** 指定震源深度和震中距采样，
        **-e** 同时计算空间导数。

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

    .. group-tab:: Python

        调用 :meth:`static_greenfn() <pygrt.pymod.PyModel1D.static_greenfn>`
        一次计算全部深度和距离组合。

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

结果保存为 ``stgrn.nc``，包含多震源深度和一维震中距采样，接收深度固定为 0 km。
静态合成自动按震源深度和源台震中距进行线性插值。
目标值不必与采样节点重合，但不能超出库的范围。

应用有限断层震源
------------------------------

**使用有限断层震源时，输出只能采用 ZNE 坐标系。**

有限震源由 |sCoulomb| 格式的断层文件描述，合成时叠加各子断层中心点源的结果。
静态解只取决于最终滑动量，不需要破裂时间函数和延迟。本例使用以下断层：

.. literalinclude:: fault.inp
    :language: text


.. figure:: fault.svg
    :align: center

    震源断层的水平投影和三维展布

本例沿走向、倾向均按 1 km 剖分，在地表网格上计算 ZNE 位移及空间导数。

.. tabs::

    .. group-tab:: CLI

        :doc:`/Module/static_syn` 模块的 **-C** 指定有限震源，**+i** 设置剖分尺寸。
        **-X/-Y** 分别指定北向、东向接收网格。

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN SRC SYN
            :end-before: END SRC SYN

    .. group-tab:: Python

        通过 :meth:`static_syn() <pygrt.pymod.PyModel1D.static_syn>` 的 *src_fault*
        读取有限震源。

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN SRC SYN
            :end-before: END SRC SYN

结果保存到 ``stsyn_ff.nc``。应力、投影及库伦应力的计算与点源相同，
这里同样取等效摩擦系数 0.4，流程见 :doc:`/Tutorial/static/static_strain_stress/static_strain_stress`。

.. tabs::

    .. group-tab:: CLI

        依次使用 :doc:`/Module/static_stress` 模块、:doc:`/Module/static_sproj` 模块和
        :doc:`/Module/static_coulomb` 模块。

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN COULOMB
            :end-before: END COULOMB

    .. group-tab:: Python

        依次调用 :func:`static_stress() <pygrt.tensor.static_stress>`、
        :func:`static_sproj() <pygrt.tensor.static_sproj>` 和
        :func:`static_coulomb() <pygrt.tensor.static_coulomb>`。

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN COULOMB
            :end-before: END COULOMB

.. figure:: disp.svg
    :align: center

    颜色表示垂直位移 Z，箭头表示水平位移 E、N，黑线表示震源断层顶边

.. figure:: coulomb.svg
    :align: center

    依次为剪应力、法向应力和库伦应力变化

应用有限断层接收
------------------------------

:doc:`/Module/static_syn` 模块和 :meth:`static_syn() <pygrt.pymod.PyModel1D.static_syn>`
也支持传入有限接收断层，本质上是批量定义多个接收点。
可以在整条断层中心接收，也可剖分后在各子断层中心接收。

**上述格林函数库的接收深度固定为 0 km。若接收断层包含地下接收点，
就必须在“接收点深度”维度上扩增格林函数库，覆盖所有接收点的深度，
形成震源深度、接收深度和震中距三个维度的采样组合。**
同时仍需确保各子源到接收点的实际震中距位于库的范围内。

接收断层与震源断层使用相同的 |sCoulomb| 文件格式，但仅支持 **Kode=100**。
滑动量不参与位移计算，行末的震源时间函数和延迟会被忽略。
几何和滑动角随结果保存，供后续应力投影使用；
若滑动角未定义，须在投影阶段补充，也可用 **rake** 格式预先给定。

同一次合成也可同时使用有限震源和有限接收断层，后处理方式与点源相同。
