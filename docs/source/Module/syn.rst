:author: 朱邓达
:date: 2025-09-22

.. include:: common_OPTs.rst_


syn
==============

:简介: 指定震源机制，根据动态格林函数合成三分量位移（及其空间导数）

语法
-----------

**grt syn** （点源）
|-G|\ *grndir*
|-S|\ [**u**]\ *scale*
|-O|\ *outdir*
[ **-Ds**\ *depsrc* ]
( [ |-R|\ *dist* ] |-A|\ *azimuth* [ **-Dr**\ *deprcv* ] | |-Q|\ *file* | |-U|\ *faultparam* )
[ |-F|\ *fn/fe/fz* | |-M|\ *strike/dip[/rake]* | |-T|\ *Mxx/Mxy/Mxz/Myy/Myz/Mzz* ]
[ |-D|\ *tftype*\ [/*tfparams*][**+d**\ *delay*] ]
[ |-I|\ *odr* ] [ |-J|\ *odr* ] [ |-N| ] [ **-e** ] [ **-s** ]
[ **-i**\ *0|1* ] [ |-P|\ *nthreads* ] [ **-h** ]

**grt syn** （有限断层）
|-G|\ *grndir*
|-C|\ *faultparam*
|-O|\ *outdir*
( [ |-R|\ *dist* ] |-A|\ *azimuth* [ **-Dr**\ *deprcv* ] | |-Q|\ *file* | |-U|\ *faultparam* )
[ |-D|\ *tftype*\ [/*tfparams*][**+d**\ *delay*] ]
[ |-I|\ *odr* ] [ |-J|\ *odr* ] [ **-e** ] [ **-s** ]
[ **-i**\ *0|1* ] [ |-P|\ *nthreads* ] [ **-h** ]


描述
--------

调用 :doc:`syn` 时，需要确定以下五类信息：

#. **格林函数输入**：|-G| 指定 :doc:`greenfn` 输出的库根目录或单个节点子目录
#. **震源位置**：点源位于水平原点，用 **-Ds** 指定深度；有限震源用 |-C| 读取断层文件
#. **台站位置**：用 **-R/-A/-Dr** 指定一个接收点，或用 |-Q| 读取任意点列表，
   或用 |-U| 读取有限接收断层，三种方式选其一
#. **震源机制与强度**：点源用 |-S| 设置强度，机制可选 |-M|、|-F| 或 |-T|，默认爆炸源；
   有限震源由 |-C| 提供位置、机制和强度，不能再设置 **-Ds/-S/-M/-F/-T**
#. **输出**：|-O| 指定 SAC 输出目录，|-N| 选择 ZNE 分量，**-e** 增加位移空间导数；
   有限震源始终输出 ZNE

输入库根目录时，**-Ds/-Dr/-R** 在库中对应维度只有一个值时可以省略。
输入单个节点子目录时，源深、台深和震中距已经确定，只需给定 |-A| 和震源机制。
有限震源或使用接收文件时须输入库根目录，所有坐标共用一个水平原点，
源台深度及源台间震中距须在库的范围内。默认使用线性插值，可用 **-i0** 改为最近邻查询。

默认合成脉冲型位移，单位为 cm；可用 |-D| 指定震源时间函数。
Z 垂直向上，R 径向向外，T 沿 R 顺时针旋转 90°，N、E 分别为北向、东向。

.. include:: explain_dynamic_output.rst_



必选选项
----------

.. _-G:

**-G**\ *grndir*
    :doc:`greenfn` 输出的库根目录或单个节点子目录。
    库根目录须包含一个模型的完整深度、距离组合及原始模型文件，
    各节点的采样点数、采样间隔和虚频率须一致，详见 :doc:`/Tutorial/dynamic/dynlib`。
    使用 |-C| 时，库根目录必须包含至少两个节点，单节点库不支持有限断层震源。
    计算位移空间导数时，所用基本源型的 **z/r** 导数文件也必须存在。

.. include:: explain_-S.rst_

.. include:: explain_-O.rst_


可选选项
--------

.. _syn-dynamic-depth:

.. _-Ds:

**-Ds**\ *depsrc*
    点源的震源深度 (km)。库中有多个震源深度时必须设置。
    输入节点子目录或使用 |-C| 时不能设置。

.. _-Dr:

**-Dr**\ *deprcv*
    单个接收点的深度 (km)。库中有多个接收深度时必须设置。
    输入节点子目录或使用 |-Q|、|-U| 时不能设置。

.. _-R:

**-R**\ *dist*
    接收点相对水平原点的距离 (km)。库中有多个震中距时必须设置。
    输入节点子目录或使用 |-Q|、|-U| 时不能设置。

.. _-A:

**-A**\ *azimuth*
    接收点相对水平原点的方位角，单位为 °，北向为 0°。
    单个接收点模式必须设置，不能与 |-Q|、|-U| 同时使用。
    范围为 [0, 360]°；源台震中距为零时，计算使用 0° 方位。

.. include:: explain_-Q.rst_

.. include:: explain_-Ufault.rst_

.. include:: explain_-Cfault.rst_

.. include:: explain_src.rst_

.. include:: explain_-Dtfunc.rst_

.. include:: explain_-IJ.rst_

.. include:: explain_rot2ZNE.rst_

.. include:: explain_-esyn.rst_

.. include:: explain_-silent.rst_

**-i**\ *0|1*
    库根目录模式下的查询方式：1 为线性插值（默认），0 为最近邻。
    输入单个节点子目录时不生效。

.. include:: explain_-P.rst_

.. include:: explain_-h.rst_


示例
-------

详见教程：

+ :doc:`/Tutorial/dynamic/syn`
