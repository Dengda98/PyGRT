:author: 朱邓达
:date: 2026-08-17

创建动态面波格林函数库
========================

动态面波解也是动态格林函数的一种。它先由 :doc:`/Module/eigenv` 计算频散结果，
再由 :doc:`/Module/modsum` 使用模态叠加法计算格林函数。由于 :doc:`/Module/modsum` 的输出目录格式与
:doc:`/Module/greenfn` **完全相同**，因此建库和合成时的深度、震中距选择方式也相同，这里只进行简单介绍。

快速上手
---------

先用 :doc:`/Module/eigenv` 模块计算 Rayleigh 波和 Love 波的频散结果，再分别调用 :doc:`/Module/modsum`。两次调用可以写入同一个输出根目录：
Rayleigh 波提供 Z、R 分量，Love 波提供 T 分量，最后得到完整的面波格林函数。

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run_library/run.sh
            :language: bash
            :start-after: BEGIN LIBRARY
            :end-before: END LIBRARY

        ``-N0`` 表示只计算基阶面波；需要更多阶数时可改为 ``-N`` 或指定阶数范围。

    .. group-tab:: Python

        .. literalinclude:: run_library/run.py
            :language: python
            :start-after: BEGIN LIBRARY
            :end-before: END LIBRARY

这里计算了 2、4 km 两个震源深度，0、2 km 两个台站深度，以及 80、100、120 km 三个震中距。
程序会遍历所有组合，目录例如 ``GRN/milrow_4_2_100/``。

.. warning::

   :doc:`/Module/eigenv` 计算频散时必须用 ``-Ff1/f2/df`` 生成等间隔频率，不能使用周期形式的 ``+p``。
   这是 :doc:`/Module/modsum` 进行逆傅里叶变换的前提。

从面波格林函数库合成
----------------------

在相同的格林函数库目录结构下，面波库与动态全波库的合成阶段完全一致，
具体用法见 :doc:`/Tutorial/dynamic/dynlib`。

若要计算空间导数及应力，建库时两个波型的 :doc:`/Module/modsum` 都须使用 **-e** 或 *calc_upar=True*。
