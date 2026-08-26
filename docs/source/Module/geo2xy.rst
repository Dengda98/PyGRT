:author: 朱邓达
:date: 2026-08-26

.. include:: common_OPTs.rst_


geo2xy
========================

:简介: 将纬度/经度坐标转换为局部 north/east 坐标


语法
-----------

**grt geo2xy**
|-G|\ *ingrid* | **-Q**\ *infile*
|-O|\ *outfile*
|-C|\ *lat0/lon0*
[ **-h** ]


描述
--------

**geo2xy** 模块使用参考点处的局部切平面近似，将地理坐标转换为局部笛卡尔坐标。
输入可以是包含地理坐标的 NetCDF 文件，也可以是坐标文本文件，输出保存到新的文件中。
输入文件不会被修改。

输入地理坐标 **lat** 和 **lon** 的单位为 degree。
输出局部坐标 **north** 和 **east** 的单位为 km。
**-C** 指定参考点的纬度和经度，格式为 *lat0/lon0*，单位为 degree。
参考点会作为新的 NetCDF 全局属性 **lat0** 和 **lon0** 写入输出文件。

逆变换公式为：

.. math::

    \mathrm{north} = (\mathrm{lat}-\mathrm{lat}_0)R\pi/180

    \mathrm{east} = \operatorname{wrap}(\mathrm{lon}-\mathrm{lon}_0)R\cos(\mathrm{lat}_0)\pi/180

其中 :math:`R=6371` km，:math:`\operatorname{wrap}` 将经度差规范化到 ``[-180, 180)``。
因此参考点在经度边界附近时，跨越 ``-180/180`` 的经度仍会得到正确的局部 east 坐标。
极点附近的坐标转换不在本模块的适用范围内。

对于 **grid** 布局，模块将 **lat** 维度和坐标变量重命名为 **north**，
将 **lon** 维度和坐标变量重命名为 **east**。
对于普通 **points** 布局和有限接收断层布局，模块保留 **point** 和 **nfault** 维度，
仅将一维坐标变量 **lat(point)** 和 **lon(point)** 重命名为 **north(point)** 和 **east(point)**。
其他变量、维度、属性和数据都会复制到输出文件中。

使用 **-Q** 时，输入文本的每个数据行前两列分别视为 **lat** 和 **lon**，
后续列按原字符串保留。
输出文件最前面会增加 ``# lat0 lon0`` 题头，用于记录参考点的纬度和经度。

必选选项
----------

.. _-G:

**-G**\ *ingrid*
    输入包含地理坐标的 NetCDF 文件。
    文件应包含 **layout** 全局属性，以及地理坐标变量 **lat** 和 **lon**。
    **-G** 与 **-Q** 不能同时使用。

.. _-Q:

**-Q**\ *infile*
    输入坐标文本文件。
    模块只读取每个数据行的前两列，分别作为 **lat** 和 **lon**。
    **-G** 与 **-Q** 不能同时使用。

.. _-O:

**-O**\ *outfile*
    输出 NetCDF 或文本文件路径。
    输出文件必须与输入文件不同。

.. _-C:

**-C**\ *lat0/lon0*
    局部坐标原点的纬度和经度，单位为 degree。
    纬度范围为 ``(-90, 90)``，经度范围为 ``[-180, 180]``。


可选选项
----------

.. include:: explain_-h.rst_
