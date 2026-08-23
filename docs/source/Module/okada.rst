:author: 朱邓达
:date: 2026-08-17

.. include:: common_OPTs.rst_


okada
==================

:简介: 使用 Okada 均匀半空间解析解计算静态位移及其空间导数

Okada 解是均匀弹性半空间中位错源产生静态变形的闭合解析解。
本模块的实现参考 NIED 的 `DC3D0/DC3D 程序说明 <https://www.bosai.go.jp/e/dc3d.html>`__，
支持埋藏点源和 Coulomb 格式有限矩形断层。关于解的背景和 PyGRT 的计算流程，
详见 :doc:`/Okada_solution/okada`。


语法
-----------

**grt okada** （点源）
|-I|\ *vp/vs/rho*
|-S|\ [**u**]\ *scale*
**-Ds**\ *depsrc*
|-O|\ *outgrid*
[ |-M|\ *strike/dip*\ [/\ *rake*] ]
[ **-Dr**\ *deprcv* ]
[ [ |-X|\ *x1/x2/dx* ] [ |-Y|\ *y1/y2/dy* ] | [ **-Q**\ *file* ] | [ **-R**\ *faultparam* ] ]
[ |-N| ] [ **-e** ] [ **-s** ] [ **-h** ]

**grt okada** （有限断层）
|-I|\ *vp/vs/rho*
|-C|\ *faultparam*
|-O|\ *outgrid*
[ **-Dr**\ *deprcv* ]
[ [ |-X|\ *x1/x2/dx* ] [ |-Y|\ *y1/y2/dy* ] | [ **-Q**\ *file* ] | [ **-R**\ *faultparam* ] ]
[ **-e** ] [ **-s** ] [ **-h** ]


描述
--------

**okada** 模块直接计算均匀弹性半空间中的静态位移。点源模式对应 Okada 的 DC3D0 解，
有限断层模式根据 Coulomb 的 ``Kode`` 解释每条记录，矩形断层使用 DC3D 解，
点源使用 DC3D0 解。

必选选项
----------

.. _-I:

**-I**\ *vp/vs/rho*
    均匀半空间参数。*vp*、*vs* 的单位为 km/s，*rho* 的单位为 g/cm\ :sup:`3`。

.. include:: explain_-S.rst_

点源模式必须设置 **-S**，有限断层模式不能使用该选项。

**-Ds**\ *depsrc*
    设置点源深度 (km)，要求非负。有限断层模式不能使用该选项，有限断层的深度从
    **-C** 文件中的顶深和底深确定。

**-Dr**\ *deprcv*
    设置规则网格接收点深度 (km)，要求非负。使用 **-Q** 时，接收点深度从文件读取，
    使用 **-R** 时接收点深度从有限断层几何读取，因此这两种模式均不能设置 **-Dr**。

.. include:: explain_-Cfault.rst_

.. include:: explain_-Ogrid.rst_

点源模式必须设置 **-S** 和 **-Ds**，有限断层模式必须设置 **-C**；两种模式不能同时使用。


可选选项
--------

.. _-M:

**-M**\ *strike/dip*\ [/\ *rake*]
    设置点源震源机制，角度单位为 °。未设置 **-M** 时为爆炸源；设置 *strike/dip* 时为张裂源；
    设置 *strike/dip/rake* 时为双力偶源。有限断层模式不能使用该选项。

.. include:: explain_-XYgrid.rst_

.. include:: explain_-Q.rst_

.. include:: explain_-Rfault.rst_

.. include:: explain_rot2ZNE.rst_

.. include:: explain_-esyn.rst_

.. include:: explain_-silent.rst_

.. include:: explain_-h.rst_


示例
-------

详见

+ :doc:`/Okada_solution/okada`
+ :doc:`/Gallery/ex18/ex18`
