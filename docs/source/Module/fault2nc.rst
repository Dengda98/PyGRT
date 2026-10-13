:author: 朱邓达
:date: 2026-10-10

.. include:: common_OPTs.rst_


fault2nc
========================

:简介: 将 |sCoulomb| 格式的有限断层文件中规则排列的剪切子断层转换为断层面 NetCDF 网格


语法
-----------

**grt fault2nc**
|-F|\ *fault*
|-O|\ *outfile*
[ **-h** ]


描述
--------

:doc:`fault2nc` 只接受满足以下条件的 |sCoulomb| 格式的有限断层文件：

* 所有记录均为 **Kode=100**，支持右旋/逆冲分量格式和 rake/net slip 格式
* 所有矩形块共面，具有相同的走向、倾角和尺寸，可拼成一个完整的大矩形
* 中心点沿走向、沿倾向分别等距排列，无缺块或重复块，输入行可乱序

输出网格的 **length**、**width** 维度分别对应沿走向、沿倾向方向，
二维变量按 **(width, length)** 存储，沿走向变化最快。
同名坐标变量和 **depth** 变量记录子断层中心位置，单位为 km。
全局属性 **strike/dip** 记录断层面走向和倾角（degree），
**top/bot** 记录整个断层面的顶深和底深（km），
**dL/dW** 记录子断层沿走向、沿倾向的尺寸（km）。
顶深取所有子断层顶深的最小值，底深取底深的最大值；

可以用 ``ncdump -h`` 命令查看输出的 |NetCDF| 网格文件，结果类似如下：

::

    netcdf delay_plane {
    dimensions:
            width = 20 ;
            length = 40 ;
    variables:
            double width(width) ;
                    width:units = "km" ;
            double length(length) ;
                    length:units = "km" ;
            double north(width, length) ;
                    north:units = "km" ;
            double east(width, length) ;
                    east:units = "km" ;
            double depth(width, length) ;
                    depth:units = "km" ;
            double rake(width, length) ;
                    rake:units = "degree" ;
            double slip(width, length) ;
                    slip:units = "m" ;
            double delay(width, length) ;
                    delay:units = "s" ;

    // global attributes:
                    :strike = 90. ;
                    :dip = 90. ;
                    :top = 1. ;
                    :bot = 21. ;
                    :dL = 5. ;
                    :dW = 1. ;
    }


必选选项
----------

.. _-F:

**-F**\ *fault*
    输入 |sCoulomb| 格式的有限断层文件。
    校验失败时直接报错，不生成输出。

.. _-O:

**-O**\ *outfile*
    输出 |NetCDF| 文件，须与输入文件不同。


示例
----------

详见 :doc:`faultdiv` 模块的示例。
