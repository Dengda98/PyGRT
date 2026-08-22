:author: 朱邓达
:date: 2026-08-17

创建静态格林函数库
===================

静态格林函数的建库思路与动态情况相同：固定模型后，把多个震源深度、台站深度和震中距交给
:doc:`/Module/static_greenfn`，由模块内部遍历全部组合计算，不需要用户手动编写循环。
区别在于，所有结果会被写入一个多维 |NetCDF| 文件，维度为
``depsrc × deprcv × north × east``，而不是生成多个文件。
保留 ``north`` 和 ``east`` 两个水平维度只是为了兼容，对应的就是震中距。


快速上手
---------

我们先简单计算这样的一个格林函数库作为示例，保存到文件 ``stgrn.nc``：

+ 震源深度：2、4 km
+ 台站深度：0、2 km
+ 震中距：0、5、10、15 km 


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

建立格林函数的目的就是为了复用。只要模型不变，以及震源深度、台站深度、震中距的间距和范围适用你的研究，就可以多次反复使用。

例如，以下我们使用刚刚计算的格林函数库，合成一个点源在二维接收平面上的结果。
与前面不同的是， **对于点源，此时的合成阶段需要明确给出震源深度，对于台站也需要给出其接收深度**。
这是因为之前章节中计算的都是单一震源深度和单一台站深度的格林函数，
如果格林函数中某一维度仅有一个，则在合成阶段可以省略指定对应参数。

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

合成位移以及位移偏导数之后，此时输出的 nc 文件格式不变，可用之前介绍的方式衔接到计算应力等步骤，
这里不再重复。

静态合成会进行插值
------------------------
仔细观察会发现，以上合成过程中， **我们指定的震源/台站深度以及震中距并没有准确的包含在静态格林函数库中**，
这正是静态库与动态库的主要使用差异。合成静态结果时，程序会进行：

1. 对包围每个网格点目标震中距的采样点，分别结合震源机制完成静态位移合成，再按目标震中距加权线性组合
2. 如果目标源深度或台站深度位于两个采样深度之间，则对各个深度组合分别合成，再按深度加权线性组合

插值作用于已经完成震源合成的位移/位移偏导数结果，而不是直接对格林函数数组插值。
目标深度以及二维网格各点的震中距必须落在库的范围内。
