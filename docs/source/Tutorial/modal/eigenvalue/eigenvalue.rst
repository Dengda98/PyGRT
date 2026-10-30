:author: 朱邓达
:date: 2026-05-20

计算频散曲线
====================

.. admonition:: 本页文档示例下载
    :class: tip

    :download:`eigenvalue.tar.gz <eigenvalue.tar.gz>` |
    解压后进入示例目录，运行 ``bash run.sh`` 即可复现本页示例。

频散曲线上每个点对应的相速度被称为面波本征值。具体地，这里我们使用自适应搜根方法 (|ars2026|) 计算面波相速度频散曲线，
命令行中使用 :doc:`/Module/eigenv` 模块，Python 中对应 :meth:`PyModel1D.eigenv() <pygrt.pymod.PyModel1D.eigenv>`。

以下示例使用如下层状模型，文件命名为 :file:`mod2` :

.. literalinclude:: mod2
    :language: text


示例程序
-----------

计算 :file:`mod2` 模型中的 Rayleigh 波和 Love 波的所有阶频散曲线，频率 5 Hz 以下，间隔 0.01 Hz。

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN DISPERSION
            :end-before: END DISPERSION

        如果仅想要快速浏览，可使用 :doc:`/Module/disp2asc` 模块将 |NetCDF|
        格式的频散转为文本格式，方便简单阅读结果以及使用 |GMT| 绘制。

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN DISPERSION
            :end-before: END DISPERSION

        可使用 :func:`read_dispersion() <pygrt.io.read_dispersion>` 读取频散结果。


计算结果分别导出为 |NetCDF| 经典格式。


.. figure:: dispersion_py.svg

绘制久期函数
---------------------

计算本征值的过程可以转变为搜索久期函数零点的过程 (|yao2026|)。
如果感兴趣，可以使用 :doc:`/Module/eigenv` 模块导出搜根过程中对应久期函数幅值。
例如，以上述程序为例，导出频率 1 Hz 处 Rayleigh 波和 Love 波的久期函数，相速度范围在 [3.0, 5.5] km/s :

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN SECFUNC
            :end-before: END SECFUNC

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN SECFUNC
            :end-before: END SECFUNC

        可使用 :func:`read_secfunc() <pygrt.io.read_secfunc>` 读取久期函数结果。


.. figure:: secfunc.svg
