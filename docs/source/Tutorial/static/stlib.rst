:author: 朱邓达
:date: 2026-08-17

创建静态格林函数库
===================

静态格林函数的建库思路与动态情况相同：固定模型后，把多个震源深度、台站深度和震中距交给
:doc:`/Module/static_greenfn`，由模块内部遍历全部组合计算，不需要用户手动编写循环。
区别在于，所有结果会被写入一个四维 |NetCDF| 文件，维度为
``depsrc × deprcv × north × east``，而不是生成多个文件。
保留 ``north`` 和 ``east`` 两个水平维度是为了兼容。


快速上手
---------

下面的命令把 2、4 km 两个震源深度，0、2 km 两个台站深度，以及 0、5、10、15 km 四个震中距保存到 ``stgrn.nc``：

.. tabs::

    .. group-tab:: C

        .. literalinclude:: run_library/run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

    .. group-tab:: Python

        .. literalinclude:: run_library/run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

        Python 接口中的深度和距离序列也必须严格递增。结果可以用
        :func:`pygrt.utils.read_static_nc` 读回，查看 ``dimensions``、``variables`` 和 ``attributes``。

静态合成会进行插值
------------------------

静态合成可以指定格林函数库中没有直接采样的目标深度和二维接收网格。
例如下面选择源深度 3 km、台站深度 1 km，并指定一个新的二维接收网格：

.. tabs::

    .. group-tab:: C

        .. literalinclude:: run_library/run.sh
            :language: bash
            :start-after: BEGIN SYN
            :end-before: END SYN

    .. group-tab:: Python

        .. literalinclude:: run_library/run.py
            :language: python
            :start-after: BEGIN SYN
            :end-before: END SYN

``3 km`` 位于 2、4 km 之间，``1 km`` 位于 0、2 km 之间，
二维网格中各点的震中距也不一定是库中的直接采样值。
这正是静态库与动态库的主要使用差异：

1. 对包围每个网格点目标震中距的采样点，分别结合震源机制完成静态位移合成，再按目标震中距加权线性组合
2. 如果目标源深度或台站深度位于两个采样深度之间，则对各个深度组合分别合成，再按深度加权线性组合

插值作用于已经完成震源合成的位移/位移偏导数结果，而不是直接对格林函数数组插值。
目标深度以及二维网格各点的震中距必须落在库的范围内。

合成时如何选择深度
------------------

:doc:`/Module/static_syn` 模块中的 **-G** 始终指向单个的静态格林函数文件。
点源合成时，多源深度库需要 **-Ds**，多台站深度库需要 **-Dr**；
只有一个对应深度时可以省略，但显式设置的值也必须与库中的值一致。
有限断层合成不设置 **-Ds**，而是根据断层几何自动使用库中的多个震源深度。

更完整的点源、任意接收点和有限断层用法见 :doc:`/Module/static_syn` 模块，Python 接口见
:meth:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` 和
:meth:`compute_static_syn() <pygrt.pymod.PyModel1D.compute_static_syn>`。
