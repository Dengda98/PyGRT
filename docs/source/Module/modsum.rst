:author: 朱邓达
:date: 2026-05-20

.. include:: common_OPTs.rst_


modsum
=================

:简介: 基于面波频散结果，使用模态叠加法计算面波格林函数

语法
-----------

**grt modsum** 
|-C|\ *path*
|-D|\ *depsrc/deprcv*
|-R|\ *r1/r2/dr*\|\ *r1,r2,...*\|\ *file*
|-O|\ *outdir*
[ |-F|\ *f1/f2* ]
[ |-N|\ [*n1*][/\ *n2*][/\ *dn*] ]
[ |-W|\ *fac* ]
[ |-E|\ [**p**]\ *t0*\ [/*v0*] ]
[ |-P|\ *nthreads* ]
[ |-G|\ **e|v|h|s** ]
[ **-e** ]
[ **-s** ]
[ **-h** ]


描述
--------

基于 **eigenv** 模块计算的相速度频散结果， **以及其结果文件中保存的模型路径，** 
**modsum** 模块使用模态叠加法计算不同震源的面波格林函数，
结果的保存方式与 :doc:`greenfn` 模块所实现的完全相同。
**modsum** 模块要求 **eigenv** 计算频散时 **-F** 必须设置等距的频率而非周期 
（即 **eigenv -F** 中不可使用 **+p**）。

必选选项
----------

.. include:: explain_-Ceigv.rst_

.. include:: explain_-D.rst_

.. include:: explain_-R.rst_

.. include:: explain_-O.rst_


可选选项
--------

.. _-F:

**-F**\ *f1/f2* 
    从频散结果中挑选频段 *[f1,f2]* 进行处理（默认全部）。

.. include:: explain_-Nmodes.rst_

.. _-W:

**-W**\ *fac*
    频率域插值倍数[1]。即在做逆傅里叶变换时，在频域上最高频后补零，
    相当于 *nt* ← *nt* \* *fac* , *dt* ← *dt* / *fac* ，使计算的波形更平滑。
    该选项与 :doc:`greenfn` **-N** 中 **+n**\ *fac* 的含义一致。

.. include:: explain_-E.rst_

.. include:: explain_-P.rst_

.. include:: explain_-G.rst_

.. include:: explain_-egrn.rst_

.. include:: explain_-s.rst_

.. include:: explain_-h.rst_


示例
-------


