:author: 朱邓达
:date: 2025-09-22
:updated_date: 2025-10-07

.. include:: common_OPTs.rst_


travt
=================

:简介: 计算半空间分层介质中的初至 P 波、S波走时

语法
-----------

**grt travt**
|-M|\ *model*
|-D|\ *depsrc/deprcv*
|-R|\ *file*\|\ *r1,r2,...*


描述
--------

**travt** 模块计算初至波的方法本质是穷举法。
在确定震源和台站相对位置的情况下，讨论“滑行”所在层位，在多种可能中比较得到走时最少的路径。
结果将打印到标准输出。


必选选项
----------

.. include:: explain_-M.rst_

.. include:: explain_-D.rst_

.. include:: explain_-R.rst_




