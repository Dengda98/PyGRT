:author: 朱邓达
:date: 2026-08-17

构建面波格林函数库与应用有限断层
========================================

.. admonition:: 本页文档示例下载
    :class: tip

    :download:`modlib.tar.gz <modlib.tar.gz>` |
    解压后进入示例目录，运行 ``bash run.sh`` 即可复现本页示例。

动态面波解也是动态格林函数的一种。它先由 :doc:`/Module/eigenv` 模块计算频散结果，
再由 :doc:`/Module/modsum` 模块使用模态叠加法计算格林函数。
:doc:`/Module/modsum` 模块与 :doc:`/Module/greenfn` 模块的输出目录格式 **完全相同**，
因此建库和合成时的深度、震中距选择方式也相同，这里只进行简单介绍。

构建面波格林函数库
------------------------

先用 :doc:`/Module/eigenv` 模块计算 Rayleigh 波和 Love 波的频散结果，再分别调用 :doc:`/Module/modsum` 模块。两次调用可以写入同一个输出根目录：
Rayleigh 波提供 Z、R 分量，Love 波提供 T 分量，最后得到完整的面波格林函数。

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN LIBRARY
            :end-before: END LIBRARY

        ``-N0`` 表示只计算基阶面波；需要更多阶数时可改为 ``-N`` 或指定阶数范围。

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN LIBRARY
            :end-before: END LIBRARY

这里计算了 2、4 km 两个震源深度，0、2 km 两个台站深度，以及 80、100、120 km 三个震中距。
程序会遍历所有组合，目录例如 ``GRN/milrow_4_2_100/``。

.. warning::

   :doc:`/Module/eigenv` 模块计算频散时必须用 ``-Ff1/f2/df`` 生成等间隔频率，不能使用周期形式的 ``+p``。
   这是 :doc:`/Module/modsum` 模块进行逆傅里叶变换的前提。

应用有限断层
----------------------

由于格林函数库的路径排布一致，因此合成的相关用法都和 :doc:`/Tutorial/dynamic/dynlib/dynlib` 中完全一致，
这里不再重复。若要计算空间导数及应力，两个波型建库时均须在 :doc:`/Module/modsum` 模块中使用 **-e**，
或在 Python 接口中设置 *calc_upar=True*。
