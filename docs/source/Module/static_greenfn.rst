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
|-O|\ *outdir*
[ |-X|\ *x1/x2/dx* ]
[ |-Y|\ *y1/y2/dy* ]
[ |-R|\ *r1/r2/dr*\|\ *r1,r2,...*\|\ *file* ]
[ |-B|\ **f|F|r|R|h|H** ]
[ |-L|\ *length*\ [**+l**\ *Flength*][**+a**\ *Ftol*][**+o**\ *offset*] ]
[ |-C|\ **d|p|n** ]
[ |-K|\ [**+k**\ *k0*][**+e**\ *keps*] ]
[ |-S| ]
[ **-e** ]
[ **-h** ]


描述
--------

**static greenfn** 模块计算多个震中距对应的静态格林函数。
但考虑到静态解的普遍应用场景，模块在计算静态格林函数时有两种传入震中距的方式：

+ 设置 |-X| 和 |-Y| 来指定二维 XY 网格。 
  坐标 (0,0) 为源点的水平投影，每个节点的震中距为 :math:`r_{ij} = \sqrt{x_i^2 + y_j^2}` 。
  实际计算中会对震中距自动去重以减少计算量。
+ 设置 |-R| 来指定多个震中距。为了兼容性，在计算和结果保存上等同于指定从原点 
  (X=0.0) 沿东向 (Y/km) 指定对应的震中距。

结果使用 NetCDF 网格保存，其中变量命名方法与 :doc:`greenfn` 模块一致。


必选选项
----------

.. include:: explain_-M.rst_

.. include:: explain_-D.rst_

.. include:: explain_-Ogrid.rst_


可选选项
--------

.. include:: explain_-XYgrid.rst_

.. include:: explain_-R.rst_

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

