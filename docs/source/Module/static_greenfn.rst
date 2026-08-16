:author: 朱邓达
:date: 2025-09-22

.. include:: common_OPTs.rst_


static_greenfn
=======================

:简介: 使用广义反射透射系数矩阵法计算静态格林函数

语法
-----------

**grt static greenfn** 
|-M|\ *model*
( |-D|\ *depsrc/deprcv* | **-Ds**\ *source* **-Dr**\ *receiver* )
|-O|\ *outgrid*
[ |-X|\ *x1/x2/dx* ]
[ |-Y|\ *y1/y2/dy* ]
[ |-R|\ *r1/r2/dr*\|\ *r1,r2,...*\|\ *file* ]
[ |-B|\ **f|F|r|R|h|H** ]
[ |-L|\ *length*\ [**+l**\ *Flength*][**+a**\ *Ftol*][**+o**\ *offset*] ]
[ |-C|\ **d|p|n** ]
[ |-K|\ [**+k**\ *k0*][**+f**][**+e**\ *keps*] ]
[ |-S| ]
[ **-e** ]
[ **-h** ]


描述
--------

**static greenfn** 模块计算静态格林函数，并将所有结果保存到一个 |NetCDF| 文件中。
文件的四个维度依次为震源深度、接收深度、北向坐标和东向坐标，
即 ``depsrc × deprcv × north × east``。保留二维水平坐标是为了兼容
|-R| 和 |-X|/|-Y| 两种输入方式。因此，一个文件可以保存多个震源深度和接收深度，
供 :doc:`static_syn` 在合成时对不同深度分别完成合成，再按目标深度加权组合结果。

建立格林函数库时，通常使用 |-R| 直接指定一维震中距列表，这是本模块的主要用法：

+ 设置 |-R| 后，程序将震中距保存为 ``north=0``、``east=R`` 的二维坐标，
  便于后续 :doc:`static_syn` 在各震中距采样点上分别完成合成，再按目标震中距加权组合合成结果。

此外，也可以使用 |-X| 和 |-Y| 指定二维网格。这种方式主要用于准确度测试，
或需要让后续合成直接使用该网格、从而避免按震中距对多个合成结果进行加权组合的场景：

+ |-X| 和 |-Y| 指定二维北向/东向网格。
  每个节点的震中距为 :math:`r_{ij} = \sqrt{north_i^2 + east_j^2}` 。
  实际计算中会对震中距自动去重以减少计算量。

静态格林函数本身只依赖震中距，而不依赖水平坐标的方向。使用 |-X|/|-Y| 时，
二维坐标用于保存测试网格，并可在后续合成时恢复接收点的方位角；使用 |-R| 时，
二维坐标只是 NC 文件兼容统一格式的一种保存方式。

震源和接收深度有两种设置方式：

+ **-D**\ *depsrc/deprcv* 设置单个震源深度和接收深度，适合兼容旧命令
+ **-Ds**\ *source* 和 **-Dr**\ *receiver* 分别设置深度列表，二者必须同时使用，
  列表语法与 |-R| 相同。多深度结果可用于点源合成时的深度加权组合和有限断层合成


必选选项
----------

.. include:: explain_-M.rst_

.. include:: explain_-D.rst_

.. note::

    **-D** 与 **-Ds/-Dr** 二选一。使用 **-Ds** 时必须同时设置 **-Dr**。

.. _-Ds:

**-Ds**\ *source*
    震源深度列表 (km)。支持逗号分隔列表、*z1/z2/dz* 等距范围和每行一个数值的文件。
    必须与 **-Dr** 一起使用，深度必须非负；程序会自动排序并合并近似重复值。

.. _-Dr:

**-Dr**\ *receiver*
    接收深度列表 (km)。语法和限制与 **-Ds** 相同，必须与 **-Ds** 一起使用。

.. include:: explain_-Ogrid.rst_


可选选项
--------

.. include:: explain_-XYgrid.rst_

.. include:: explain_-R.rst_

.. include:: explain_-Bbound.rst_

.. include:: explain_-L.rst_

.. _-K:

**-K**\ [**+k**\ *k0*][**+f**][**+e**\ *keps*]
    控制波数积分搜索区间的上界
    :math:`k_{\text{max,ref}} = k_0 \cdot \dfrac{\pi}{\Delta h}`，其中
    :math:`\Delta h = \max(|z_s-z_r|, 0.1)` km。

    程序在 :math:`[\Delta k, k_{\text{max,ref}}]` 内基于核函数振幅搜索实际积分上限 :math:`k_{\text{max}}` 。
    若搜索达到 :math:`k_{\text{max,ref}}` 仍未收敛，或震源与场点完全同深度时，
    默认模式下将自动启用 DCM 。

    + **+k**\ *k0* - 零频项系数 [50.0]，
      其中深度差 :math:`\Delta h = \max(|z_s - z_r|, 0.1)` 。
    + **+f** - 直接使用 :math:`k_{\text{max,ref}}` 作为积分上限，
      不进行振幅搜索。
    + **+e**\ *keps* - 用于判断提前结束波数积分的收敛精度[0.0, 默认不使用]，
      详见 Yao and Harkrider (1983) 和 :doc:`/Advanced/k_integ/kmax` 。

    当震源和接收点同深度，或自动搜索达到参考上限仍未收敛时，自动收敛模式会使用 DCM。

.. include:: explain_-Cconverg.rst_

.. _-S:

**-S**
    输出波数积分过程中的核函数文件，保存目录为 ``stgrtstats`` 。该选项仅适用于
    单个震源深度和单个接收深度；计算多深度格林函数时不会输出这些统计文件。
    关于文件格式及其读取详见 :doc:`/Advanced/integ_converg/integ_converg` 。

    .. include:: explain_-Sstats.rst_

.. include:: explain_-egrn.rst_

.. include:: explain_-h.rst_



示例
-------

详见教程：

+ :doc:`/Tutorial/static/static_gfunc`
+ :doc:`/Advanced/integ_converg/integ_converg`
