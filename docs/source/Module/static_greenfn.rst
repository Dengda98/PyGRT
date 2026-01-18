:author: 朱邓达
:date: 2025-09-22

.. include:: common_OPTs.rst_


static_greenfn
=======================

:简介: 使用广义反射透射系数矩阵法 + 离散波数法计算静态格林函数

语法
-----------

**grt static greenfn** 
|-M|\ *model*
|-D|\ *depsrc/deprcv*
|-X|\ *x1/x2/dx*
|-Y|\ *y1/y2/dy*
|-O|\ *outdir*
[ |-B|\ **f|F|r|R|h|H** ]
[ |-L|\ *length*\ [**+l**\ *Flength*][**+a**\ *Ftol*][**+o**\ *offset*] ]
[ |-C|\ **d|p|n** ]
[ |-K|\ [**+k**\ *k0*][**+e**\ *keps*] ]
[ |-S| ]
[ **-e** ]
[ **-h** ]


描述
--------

**static greenfn** 模块将震源置于原点 *(0,0)* ，计算由 |-X| （北向）和 |-Y| （东向）
定义的 2D 网格点上的静态格林函数。与动态情况的 :doc:`greenfn` 模块指定 1D 的震中距数组不同，
这里使用 2D 网格可能会在某些特殊的对称情况下增加重复计算，但会方便其它模块的数据交互。
输出结果使用 NetCDF 网格保存，其中变量命名方法与 :doc:`greenfn` 模块一致。


必选选项
----------

.. include:: explain_-M.rst_

.. include:: explain_-D.rst_

.. _-X:

**-X**\ *x1/x2/dx*
    指定北向等距网格点。 *x1,x2* 分别为最小、最大值， *dx* 为间隔。

.. _-Y:

**-Y**\ *y1/y2/dy*
    指定东向等距网格点。 *y1,y2* 分别为最小、最大值， *dy* 为间隔。

.. include:: explain_-Ogrid.rst_



可选选项
--------

.. include:: explain_-Bbound.rst_

.. include:: explain_-L.rst_

.. _-K:

**-K**\ [**+k**\ *k0*][**+e**\ *keps*]
    控制波数积分上限 :math:`k_0 \cdot \dfrac{\pi}{\Delta h}`

    + **+k**\ *k0* - 控制零频的积分上限 [5.0]，其中深度差 :math:`\Delta h = \max(|z_s - z_r|, 1.0)` 。
    + **+e**\ *keps* - 用于判断提前结束波数积分的收敛精度[0.0, 默认不使用]，详见
      Yao and Harkrider (1983) 和 :doc:`/Advanced/k_integ/kmax` 。

.. include:: explain_-Cconverg.rst_

.. _-S:

**-S**
    输出波数积分过程中的核函数文件，保存目录为 ``stgrtstats`` 。
    关于文件格式及其读取详见 :doc:`/Advanced/integ_converg/integ_converg` 。

    .. include:: explain_-Sstats.rst_

.. include:: explain_-egrn.rst_

.. include:: explain_-h.rst_



示例
-------

详见教程：

+ :doc:`/Tutorial/static/static_gfunc`
+ :doc:`/Advanced/integ_converg/integ_converg`

