.. include:: common_OPTs.rst_


syn
==============

:简介: 指定震源机制，根据动态格林函数合成三分量位移（及其空间导数）

语法
-----------

**grt syn** 
|-G|\ *grndir*
|-A|\ *azimuth*
|-S|\ [**u**]\ *scale*
|-O|\ *outdir*
[ |-F|\ *fn/fe/fz* ]
[ |-M|\ *strike/dip/rake* ]
[ |-T|\ *Mxx/Mxy/Mxz/Myy/Myz/Mzz* ]
[ |-D|\ *tftype/tfparams* ]
[ |-I|\ *odr* ]
[ |-J|\ *odr* ]
[ |-N| ]
[ **-e** ]
[ **-s** ]
[ **-h** ]


描述
--------

**syn** 模块默认合成的是脉冲型 （impulse-like）位移，单位 cm。三分量位移分别为：

+ **Z** - 垂直向上
+ **R** - 径向向外
+ **T** - 沿 **R** 方向顺时针旋转 90 °

保存 SAC 文件路径为 ``*{outdir}/{ch}.sac*``，其中 ``ch`` 为分量名。
若设置 |-N| 则输出为 ZNE 分量，其中 Z 不变，NE分别为北向和东向。

目前合成位移支持的震源类型包括爆炸源(|-S|)，单力源(|-F|)，剪切源(|-M|)和矩张量源(|-T|)，
其中矩张量源是根据爆炸源和剪切源合成的。一次只能选择一种震源。



必选选项
----------

.. _-G:
    
**-G**\ *grndir*
    保存格林函数的目录，例如 :doc:`greenfn` 模块计算的单个震中距的格林函数保存在
    ``{outdir}/{model}_{depsrc}_{deprcv}_{r}`` 。

.. _-A:

**-A**\ *azimuth*
    方位角，单位°，北向为 0 °

.. include:: explain_-S.rst_

.. include:: explain_-O.rst_



可选选项
--------

.. include:: explain_src.rst_

.. include:: explain_-Dtfunc.rst_

.. include:: explain_-IJ.rst_

.. include:: explain_rot2ZNE.rst_

.. include:: explain_-esyn.rst_

.. include:: explain_-s.rst_

.. include:: explain_-h.rst_



示例
-------

详见教程：

+ :doc:`/Tutorial/dynamic/syn`
