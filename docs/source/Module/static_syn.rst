:author: 朱邓达
:date: 2025-09-22

.. include:: common_OPTs.rst_


static_syn
==================

:简介: 指定震源机制，根据静态格林函数合成三分量位移（及其空间导数）


语法
-----------

**grt static syn** （点源）
|-G|\ *ingrid*
|-S|\ [**u**]\ *scale*
|-O|\ *outgrid*
[ **-Ds**\ *depsrc* ] [ **-Dr**\ *deprcv* ]
[ |-F|\ *fn/fe/fz* | |-M|\ *strike/dip[/rake]* | |-T|\ *Mxx/Mxy/Mxz/Myy/Myz/Mzz* ]
[ [ |-X|\ *x1/x2/dx* ] [ |-Y|\ *y1/y2/dy* ] | [ **-Q**\ *file* ] | [ **-U**\ *faultparam* ] ]
[ |-N| ]
[ |-P|\ *nthreads* ]
[ **-e** ]
[ **-s** ]
[ **-h** ]

**grt static syn** （有限断层）
|-G|\ *ingrid*
|-C|\ *faultparam*
|-O|\ *outgrid*
[ **-Dr**\ *deprcv* ]
[ [ |-X|\ *x1/x2/dx* ] [ |-Y|\ *y1/y2/dy* ] | [ **-Q**\ *file* ] | [ **-U**\ *faultparam* ] ]
[ |-P|\ *nthreads* ]
[ **-e** ] [ **-s** ] [ **-h** ]


描述
--------

调用 :doc:`static_syn` 时，需要确定以下五类信息：

#. **格林函数输入**：|-G| 指定 :doc:`static_greenfn` 生成的单个四维 |NetCDF| 格林函数库
#. **震源位置**：点源位于水平原点，用 **-Ds** 指定深度；有限震源用 |-C| 从断层文件读取各震源位置
#. **台站位置**：默认沿用库中的水平网格，也可用 **-X/-Y** 重新定义网格；
   网格上的台站共用 **-Dr** 指定的深度。另可选择 |-Q| 读取逐点坐标，
   或用 |-U| 读取有限接收断层；接收文件已经提供各点深度，不能再设置 **-X/-Y/-Dr**
#. **震源机制与强度**：点源必须设置 |-S|，再至多选择 |-M|、|-F|、|-T| 中的一种机制；
   均不选择时为爆炸源。|-M| 给出走向、倾角和滑动角时为剪切源，省略滑动角时为张裂源。
   有限震源的 |-C| 同时提供位置、机制和滑动量或矩势，因此不再设置 **-Ds/-S/-M/-F/-T**
#. **输出**：|-O| 指定 |NetCDF| 输出文件；|-N| 选择 ZNE 分量，**-e** 增加位移空间导数。
   有限震源始终输出 ZNE

点源的 **-Ds** 和网格台站的 **-Dr**，在库中对应深度维度只有一个值时可以省略。
所有坐标共用一个水平原点，源台几何应位于库的采样范围内。

库采样点之间的结果采用插值：先在相邻距离及深度节点上结合震源机制完成合成，
再按目标位置加权组合。位移单位为 cm，默认 Z 垂直向上，R 径向向外，T 沿 R 顺时针旋转 90°。

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

.. _-G:
    
**-G**\ *ingrid*
    :doc:`static_greenfn` 生成的静态格林函数库文件。必须是单个四维 |NetCDF| 文件。

.. include:: explain_-S.rst_

.. _static-syn-depth:

**-Ds**\ *depsrc*
    点源的震源深度 (km)。当输入格林函数库含多个震源深度时必须设置；
    只有一个震源深度时可以省略，也可以显式设置，但显式设置的值必须与库中深度一致；
    有限断层模式禁止设置。深度必须落在格林函数库的震源深度范围内。
    如果位于两个采样深度之间，程序会分别使用两侧深度完成合成，再按目标深度加权组合结果。

**-Dr**\ *deprcv*
    网格接收点的接收深度 (km)。当格林函数库含多个接收深度且未使用 **-Q** 或 **-U** 时必须设置；
    输入库只有一个接收深度时可以省略，也可以显式设置，但显式设置的值必须与库中深度一致。
    使用 **-Q** 时接收深度从文件读取，因此仍禁止设置 **-Dr**。
    使用 **-U** 时接收深度从有限断层几何读取，因此也禁止设置 **-Dr**。
    如果目标接收深度位于两个采样深度之间，程序会分别使用两侧深度完成合成，
    再按目标接收深度加权组合结果。

.. include:: explain_-Ogrid.rst_



可选选项
--------

.. include:: explain_-XYgrid.rst_

.. include:: explain_-Q.rst_

.. include:: explain_-Ufault.rst_

.. include:: explain_-Cfault.rst_

.. include:: explain_src.rst_

.. include:: explain_rot2ZNE.rst_

.. include:: explain_-esyn.rst_

.. include:: explain_-P.rst_

.. include:: explain_-silent.rst_

.. include:: explain_-h.rst_



示例
-------

详见教程：

+ :doc:`/Tutorial/static/static_syn`
