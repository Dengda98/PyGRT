:author: 朱邓达
:date: 2025-04-17

合成静态位移
=====================

Python中合成静态位移的主函数为 :meth:`static_syn() <pygrt.pymod.PyModel1D.static_syn>` ，
C模块为 :doc:`/Module/static_syn`。

使用上节计算的格林函数，我们可以指定 **点源的震源机制** 以及 **接收点位置** 来合成静态位移。
为方便画图，以下结果都使用ZNE分量。
若仅使用已算好的静态格林函数，构造 :class:`~pygrt.pymod.PyModel1D` 时只需指定 ``stgrn``，无需再传 ``modelpath``：

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN REUSE STGRN
    :end-before: END REUSE STGRN

二维网格示例使用 |GMT| 绘图，任意点示例使用 Python 绘图。

这里提供计算和绘图的 Python 脚本和 Shell 脚本供下载参考。
:download:`Shell Scripts <run/run.sh>` | :download:`Python Scripts <run/run.py>`


指定二维接收网格
---------------------
我们可以通过设置北向坐标和东向坐标来定义一个二维的水平网格点，然后计算这些规则网格点上的静态位移。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN GRID
            :end-before: END SYN GRID

        可以使用 **-X, -Y** 来指定二维接收网格。
        各种点源的震源参数设置用法详见 :doc:`/Module/static_syn` 模块。

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN GRID
            :end-before: END SYN GRID

        可以设置 **norths, easts** 来指定二维接收网格。
        各种点源的震源参数设置用法详见 
        :meth:`static_syn() <pygrt.pymod.PyModel1D.static_syn>`
        函数 API 说明。


.. figure:: run/syn_grid.svg
    :align: center

    黑色箭头表示水平位移，颜色表示垂直位移


指定任意点
-------------
我们也可以通过给定一个记录有任意点坐标的文件来计算这些点上的位移，
这个文件要求前三列分别为北向坐标、东向坐标以及深度，单位均为 km，
以 “#” 开头的行会被忽略。

下面只以一个走滑源 -M33/90/0 为例。
我们示意性的在震中水平面上取半径为 5 km 的圆环，
均匀生成 72 个接收点，然后计算各点的静态位移。
此时该文件内部大概看起来像这样：

.. literalinclude:: run/rcv_head
    :language: text

.. tabs::

    .. group-tab:: C

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN POINTS
            :end-before: END SYN POINTS

    .. group-tab:: Python

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN POINTS
            :end-before: END SYN POINTS

绘图脚本仅供参考（:download:`plot_points.py <run/plot_points.py>`），
实际研究中，点的分布不同（例如通常会要计算沿着某条断层线上的位移以及衍生的物理量），
合适的呈现形式也要进行修改。

.. figure:: run/syn_points.svg
    :align: center

    黑色箭头表示水平位移，颜色表示垂直位移
