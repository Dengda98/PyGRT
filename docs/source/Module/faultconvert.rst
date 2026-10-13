:author: 朱邓达
:date: 2026-10-12

.. include:: common_OPTs.rst_


faultconvert
========================

:简介: 将 USGS/WISP 的 PARAM、FSP 格式的动态有限断层转换为 |sCoulomb| 格式


语法
-----------

**grt faultconvert**
|-F|\ **[p|f]/**\ *fault*
|-C|\ *lat/lon*
|-O|\ *file*
[ |-N|\ *index* ]
[ **-h** ]


描述
--------

:doc:`faultconvert` 将 USGS/WISP 发布的 PARAM、FSP 动态有限断层转换为 |sCoulomb| 格式，
支持单段和多段断层。两种输入格式可从 USGS 上选择地震事件来下载示例，例如
`PARAM 文件 <https://earthquake.usgs.gov/product/finite-fault/us6000tkt2_1/us/1787347333870/basic_inversion.param>`_
和 `FSP 文件 <https://earthquake.usgs.gov/product/finite-fault/us6000tkt2_1/us/1787347333870/complete_inversion.fsp>`_。

转换结果保留子断层的滑动量和破裂过程，可用于 :doc:`syn` 或 :doc:`lamb` 合成。

USGS/WISP 的 PARAM 数据行保存的是每块子断层的 **首个内部点坐标**，通常偏离中心；
FSP 数据行保存的则是子断层中心坐标。PARAM 的第三个边界点和末条数据记录，
分别给出最后一块子断层的末个、首个内部点，两者经纬度、深度差值的一半就是中心偏移。
本模块将该偏移加到 PARAM 首条记录上，恢复首块子断层中心。

每个断层段只投影首块子断层的中心，再按走向、倾角、标称剖分尺寸和导出顺序重建规则矩形。
两种格式都按沿走向索引优先变化的顺序展开。


必选选项
----------

.. _-F:

**-F[p|f]/**\ *fault*
    输入断层文件，必须显式指定格式：

    + **p** - PARAM 格式文件，包含破裂速率的上升时长和下降时长
    + **f** - FSP 格式文件，破裂速率的上升时长和下降时长被合并为总时长，此时本模块仅作上升下降时长等分的处理。

.. _-C:

**-C**\ *lat/lon*
    局部坐标原点的纬度和经度，单位为 degree。
    纬度范围为 (−90, 90)，经度范围为 [−180, 180]。
    用于将地理坐标转为局部笛卡尔坐标。

.. _-O:

**-O**\ *file*
    输出转换的 |sCoulomb| 文件。


可选选项
----------

.. _-N:

**-N**\ *index*
    仅转换指定断层段，按文件中的段顺序从 1 编号，须为正整数。
    省略时，将全部断层段依次写入同一个输出文件。


示例
----------

+ :doc:`/Gallery/ex23/ex23`
