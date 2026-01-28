:author: 朱邓达
:date: 2025-04-17

合成静态位移
=================

Python中合成静态位移的主函数为 :func:`gen_syn_from_gf_*() <pygrt.utils.gen_syn_from_gf_DC>` （\*表示对不同震源）（与合成动态位移的函数共用）  ，C模块为 :doc:`/Module/static_syn`。

使用上节计算的格林函数，合成静态位移。为方便画图，以下结果都使用ZNE分量。

不同震源
-------------
以下绘图使用 `GMT <https://www.generic-mapping-tools.org/>`_  绘制。这里提供计算和绘图的 Python 脚本和 Shell 脚本供下载参考。

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

程序也支持在合成阶段指定新的 XY 网格，此时每个节点将使用最近震中距的格林函数来进行合成，
因此这是一种近似方法，程序会打印震中距之差的统计信息，以方便评估近似效果。
但近似好处是一旦静态格林函数计算好，合成阶段使用新的网格也可以复用格林函数，无需再计算。

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


