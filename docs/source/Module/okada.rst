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
**-Ds**\ *depsrc*
|-S|\ [**u**]\ *scale*
|-O|\ *outgrid*
( |-X|\ *x1/x2/dx* |-Y|\ *y1/y2/dy* **-Dr**\ *deprcv* | |-Q|\ *file* | |-U|\ *faultparam* )
[ |-M|\ *strike/dip*\ [/\ *rake*] ]
[ |-N| ] [ **-e** ] [ **-s** ] [ **-h** ]

**grt okada** （有限断层）
|-I|\ *vp/vs/rho*
|-C|\ *faultparam*
|-O|\ *outgrid*
( |-X|\ *x1/x2/dx* |-Y|\ *y1/y2/dy* **-Dr**\ *deprcv* | |-Q|\ *file* | |-U|\ *faultparam* )
[ |-N| ] [ **-e** ] [ **-s** ] [ **-h** ]


描述
--------

:doc:`okada` 模块直接计算均匀弹性半空间中的静态位移。点源模式对应 Okada 的 DC3D0 解，
有限断层模式根据 Coulomb 的 *Kode* 字段解释每条记录，矩形断层使用 DC3D 解，
点源使用 DC3D0 解。

调用 :doc:`okada` 时，需要确定以下五类信息：

#. **介质**：|-I| 指定均匀半空间参数，无需格林函数库
#. **震源位置**：点源位于水平原点，必须用 **-Ds** 指定深度；有限震源用 |-C| 读取断层文件
#. **接收位置**：用 **-X/-Y/-Dr** 指定规则网格，三个选项都需设置；
   或用 |-Q| 读取任意点列表，或用 |-U| 读取有限接收断层，三种方式选其一
#. **震源机制与强度**：点源用 |-S| 设置强度，机制由 |-M| 指定，默认爆炸源；
   有限震源由 |-C| 提供位置、机制和强度，不能再设置 **-Ds/-S/-M**
#. **输出**：|-O| 指定 NetCDF 输出文件，|-N| 选择 ZNE 分量，**-e** 增加位移空间导数

输出文件的全局 **layout** 属性区分三种接收布局：

* **grid**：规则网格，结果使用 **north/east** 二维布局
* **points**：|-Q| 指定的任意点列表，结果使用一维 **point** 维度；
  若提供接收断层形态，还会保存逐点 **strike**、**dip**、**rake** 变量
* **faults**：|-U| 指定的有限接收断层，结果使用一维 **point** 维度，
  并保存断层级形态和剖分信息。**nfault** 为断层数量，**offset** 为各断层接收点的结束索引，
  **stksize**、**dipsize** 分别为沿走向、倾向的子断层数量

**points** 和 **faults** 均保存各点的 **north**、**east** 和 **depth** 坐标。

必选选项
----------

.. _-I:

**-I**\ *vp/vs/rho*
    均匀半空间参数。*vp*、*vs* 的单位为 km/s，*rho* 的单位为 g/cm\ :sup:`3`。

.. include:: explain_-S.rst_

.. include:: explain_-Ogrid.rst_


可选选项
--------

.. _-Ds:

**-Ds**\ *depsrc*
    点源的震源深度 (km)，允许取 0。点源模式必须设置，使用 |-C| 时不能设置。

.. _-Dr:

**-Dr**\ *deprcv*
    规则网格接收点的深度 (km)，允许取 0。网格模式必须设置，使用 |-Q|、|-U| 时不能设置。

.. include:: explain_-XYgrid.rst_

.. include:: explain_-Q.rst_

.. include:: explain_-Ufault.rst_

.. include:: explain_-Cfault.rst_


.. _-M:

**-M**\ *strike/dip*\ [/\ *rake*]
    设置点源震源机制，角度单位为 °。未设置 |-M| 时为爆炸源；设置 *strike/dip* 时为张裂源；
    设置 *strike/dip/rake* 时为双力偶源。使用 |-C| 时不能设置。

.. include:: explain_rot2ZNE.rst_

.. include:: explain_-esyn.rst_

.. include:: explain_-silent.rst_

.. include:: explain_-h.rst_


示例
-------

详见

+ :doc:`/Okada_solution/okada`
+ :doc:`/Gallery/ex18/ex18`
