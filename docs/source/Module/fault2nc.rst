:author: 朱邓达
:date: 2026-10-10

.. include:: common_OPTs.rst_


fault2nc
========================

:简介: 将 |simplified-Coulomb| 格式的有限断层文件中规则排列的剪切子断层转换为断层面 NetCDF 网格


语法
-----------

**grt fault2nc**
|-F|\ *fault*
|-O|\ *outfile*
[ **-h** ]


描述
--------

:doc:`fault2nc` 只接受满足以下条件的 |simplified-Coulomb| 格式的有限断层文件：

* 所有记录均为 **Kode=100**，支持右旋/逆冲分量格式和 rake/net slip 格式
* 所有矩形块共面，具有相同的走向、倾角和尺寸，可拼成一个完整的大矩形
* 中心点沿走向、沿倾向分别等距排列，无缺块或重复块，输入行可乱序

可以用 ``ncdump -h`` 命令查看输出的 |NetCDF| 网格文件，结果类似如下：

::

    netcdf delay_plane {
    dimensions:
            dip = 20 ;
            strike = 40 ;
    variables:
            double dip(dip) ;
                    dip:units = "km" ;
            double strike(strike) ;
                    strike:units = "km" ;
            double north(dip, strike) ;
                    north:units = "km" ;
            double east(dip, strike) ;
                    east:units = "km" ;
            double depth(dip, strike) ;
                    depth:units = "km" ;
            double rake(dip, strike) ;
                    rake:units = "degree" ;
            double slip(dip, strike) ;
                    slip:units = "m" ;
            double delay(dip, strike) ;
                    delay:units = "s" ;

    // global attributes:
                    :strike = 90. ;
                    :dip = 90. ;
    }


必选选项
----------

.. _-F:

**-F**\ *fault*
    输入 |simplified-Coulomb| 格式的有限断层文件。
    校验失败时直接报错，不生成输出。

.. _-O:

**-O**\ *outfile*
    输出 |NetCDF| 文件，须与输入文件不同。


示例
----------

详见 :doc:`faultdiv` 模块的示例。
