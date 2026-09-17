:author: 朱邓达
:date: 2026-09-16

.. include:: common_OPTs.rst_


lamb
==============

:简介: 组合三类 Lamb 闭合解，直接计算均匀弹性半空间中的动态全波解

语法
-----------

**grt lamb**
|-H|\ *vp/vs/rho*
|-N|\ *nt/dt*
|-R|\ *dist*
**-Ds**\ *depsrc*
**-Dr**\ *deprcv*
|-A|\ *azimuth*
|-S|\ [**u**]\ *scale*
|-O|\ *outdir*
[ |-F|\ *fn/fe/fz* ]
[ |-M|\ *strike/dip[/rake]* ]
[ |-T|\ *Mxx/Mxy/Mxz/Myy/Myz/Mzz* ]
[ |-D|\ *tftype/tfparams* ]
[ |-E|\ [**p**]\ *t0*\ [/*v0*] ]
[ |-I|\ *odr* ]
[ |-J|\ *odr* ]
[ **-n** ]
[ **-e** ]
[ **-s** ]
[ **-h** ]


描述
--------

**lamb** 模块将第一、二、三类 Lamb 问题的广义闭合解组合成一个统一接口，
根据震源深度和接收点深度自动选择对应的计算方法：

+ 两点都在自由表面时调用第一类 Lamb 解（**lamb1**）
+ 只有一个点位于地下时调用第二类 Lamb 解（**lamb2**）
+ 两点都位于地下时调用第三类 Lamb 解（**lamb3**）

**lamb** 模块负责将 Lamb 闭合解的无量纲结果恢复为物理量，
再按照震源机制进行组合，并以与 :doc:`syn` 类似的 SAC 文件形式输出。
各个震相的到时和名称也将写入 SAC 头段变量中。

设置 **-n** 后输出为 ZNE 分量。设置 **-e** 后，同时输出位移的空间导数.

两点都在自由表面时，第一类 Lamb 解只支持单力源 **-F**，此时 **-e** 会被忽略。


必选选项
----------

.. _-H:

**-H**\ *vp/vs/rho*
    均匀半空间参数。*vp* 和 *vs* 的单位为 km/s，*rho* 的单位为 g/cm\ :sup:`3`。

.. _-N:

**-N**\ *nt/dt*
    输出采样点数 *nt* 和采样间隔 *dt*，*dt* 的单位为 s。记录从发震时刻开始。

.. _-R:

**-R**\ *dist*
    源点和接收点之间的水平震中距，单位为 km，要求为正数。

.. _-Ds:

**-Ds**\ *depsrc*
    震源深度，单位为 km，允许取 0。

.. _-Dr:

**-Dr**\ *deprcv*
    接收点深度，单位为 km，允许取 0。

.. _-A:

**-A**\ *azimuth*
    从震源指向接收点的方位角，单位为 °，北向为 0 °，取值范围为 [0, 360]。

.. _-S:

**-S**\ [**u**]\ *scale*
    震源放大系数。爆炸源、剪切源和矩张量源中 *scale* 为标量地震矩，单位为 dyne-cm；
    单力源中 *scale* 的单位为 dyne。添加 **u** 后，*scale* 被解释为震源面积乘滑动量，
    程序会再乘以震源处的剪切模量。

.. _-O:

**-O**\ *outdir*
    输出目录名，不存在时自动新建。已有同名文件会被覆盖。


可选选项
--------

.. include:: explain_src.rst_

.. include:: explain_-Dtfunc.rst_

.. include:: explain_-E.rst_

.. include:: explain_-IJ.rst_

.. _lamb-n:

**-n**
    将输出分量设置为 ZNE，默认为 ZRT。此处的选项名为小写 **-n**，
    与 :doc:`syn` 中用于同一功能的大写 **-N** 不同。

.. _lamb-e:

**-e**
    同时计算位移空间导数。输出文件名在分量名前增加导数方向前缀，
    例如 ``zN.sac`` 表示北向位移对接收点 z 坐标的偏导。

.. include:: explain_-silent.rst_

.. include:: explain_-h.rst_


示例
-------

+ :doc:`/Lamb_problem/lamb`
