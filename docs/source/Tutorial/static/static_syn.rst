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


.. image:: run/syn_ex.png
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


.. image:: run/syn_sf.png
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


.. image:: run/syn_dc.png
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


.. image:: run/syn_dc2.png
    :width: 500px
    :align: center




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


.. image:: run/syn_mt.png
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


.. image:: run/syn_mt2.png
    :width: 500px
    :align: center