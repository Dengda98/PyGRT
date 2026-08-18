:author: 朱邓达
:date: 2025-04-17

合成静态位移
=================

Python中合成静态位移的主函数为 :meth:`compute_static_syn() <pygrt.pymod.PyModel1D.compute_static_syn>` ，
C模块为 :doc:`/Module/static_syn`。

使用上节计算的格林函数，合成静态位移。为方便画图，以下结果都使用ZNE分量。
若仅使用已算好的静态格林函数，构造 :class:`~pygrt.pymod.PyModel1D` 时只需指定 ``stgrn``，无需再传 ``modelpath``：

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN REUSE STGRN
    :end-before: END REUSE STGRN

不同震源
-------------
Python 接口根据震源专用参数自动确定震源类型：不设置 ``strike``、``dip``、``rake``、
``force`` 和 ``moment_tensor`` 时为爆炸源；设置
``force`` 时为单力源；设置 ``moment_tensor`` 时为矩张量源；设置 ``strike`` 和
``dip`` 时为张裂源，若同时设置 ``rake`` 则为剪切源。一次只能设置一组震源专用参数，
不完整或混用参数会报错。有限断层使用 ``finite_fault`` 参数，与点源参数互斥。

以下绘图使用 |GMT| 绘制。这里提供计算和绘图的 Python 脚本和 Shell 脚本供下载参考。

:download:`Shell Scripts <run/run.sh>`

:download:`Python Scripts <run/run.py>`


爆炸源
~~~~~~~~~~~~~~~~~
标量矩 1e24 dyne·cm。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN EX
            :end-before: END SYN EX

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN EX
            :end-before: END SYN EX


.. figure:: run/syn_ex.svg
    :width: 500px
    :align: center


单力源
~~~~~~~~~~~~~~~~~
北向力 :math:`f_N=1`，东向力 :math:`f_E=-0.5`，垂直向下的力 :math:`f_Z=2`，单位 1e16 dyne。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN SF
            :end-before: END SYN SF

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN SF
            :end-before: END SYN SF


.. figure:: run/syn_sf.svg
    :width: 500px
    :align: center


剪切源
~~~~~~~~~~~~~~
断层走向33°，倾角50°，滑动角120°，标量矩 1e24 dyne·cm。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN DC
            :end-before: END SYN DC

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN DC
            :end-before: END SYN DC


.. figure:: run/syn_dc.svg
    :width: 500px
    :align: center

这里如果改变倾角为90°，滑动角0°，就可以看到清晰的蝴蝶状辐射花样。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN DC2
            :end-before: END SYN DC2

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN DC2
            :end-before: END SYN DC2


.. figure:: run/syn_dc2.svg
    :width: 500px
    :align: center


张裂源
~~~~~~~~~~~~~~
断层走向33°，倾角50°，标量矩 1e24 dyne·cm。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN TS
            :end-before: END SYN TS

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN TS
            :end-before: END SYN TS


.. figure:: run/syn_ts.svg
    :width: 500px
    :align: center

    沙滩球仅绘制了张裂源中的 DC+CLVD 分量


这里如果改变倾角为90°，就可以看到清晰的对称辐射花样。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN TS2
            :end-before: END SYN TS2

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN TS2
            :end-before: END SYN TS2


.. figure:: run/syn_ts2.svg
    :width: 500px
    :align: center

    沙滩球仅绘制了张裂源中的 DC+CLVD 分量


矩张量源
~~~~~~~~~~~~~~
:math:`M_{xx}=0.1, M_{xy}=-0.2, M_{xz}=1.0, M_{yy}=0.3, M_{yz}=-0.5, M_{zz}=-2.0`，单位 1e24 dyne·cm， **其中X为北向，Y为东向，Z为垂直向下**。


.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN MT
            :end-before: END SYN MT

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN MT
            :end-before: END SYN MT


.. figure:: run/syn_mt.svg
    :width: 500px
    :align: center


若指定 :math:`M_{xy}=-0.2`，其它均为零，则为纯剪切。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SYN MT2
            :end-before: END SYN MT2

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SYN MT2
            :end-before: END SYN MT2


.. figure:: run/syn_mt2.svg
    :width: 500px
    :align: center


指定新的 XY 网格
---------------------

以上示例结果中， XY 网格位置都沿用了静态格林函数计算时传入的 XY 网格。

程序也支持在合成阶段指定新的 XY 网格。对于每个新网格节点，程序先使用格林函数库中包围
该节点震中距的采样点分别完成合成，再根据该节点的震中距计算权重并组合这些合成结果；
如果还需要根据深度进行权重组合，则先对包围目标深度的各个深度组合分别完成上述合成，再根据目标深度
计算权重并组合结果。这里组合的是已经完成震源合成的结果，而不是直接对格林函数数组插值。
只要新的网格中的震中距和深度落在格林函数库范围内，静态格林函数计算好后，合成阶段就可以
复用该格林函数库，无需重新计算。

以下以一个与上面相同的走滑断层作为示例进行计算，选取了间隔稍大的网格。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN NEW XY
            :end-before: END NEW XY

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN NEW XY
            :end-before: END NEW XY

.. grid:: 2

    .. grid-item::

        .. figure:: run/syn_dc2.svg
            :width: 400px
            :align: center

    .. grid-item::

        .. figure:: run/synXY_dc2.svg
            :width: 400px
            :align: center
