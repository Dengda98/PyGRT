:author: 朱邓达
:date: 2025-04-17

计算静态格林函数
=================

Python中计算静态格林函数的主函数为 :meth:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` ，C模块为 :doc:`/Module/static_greenfn`。

建议先阅读完 :doc:`/Tutorial/dynamic/gfunc` 部分。静态情况与动态情况采取的计算方法一致，只是推导细节会有不同，详见 |yao2026p| 。

静态解模块对于传入“震中距”的方式及后续处理与动态解模块相比稍有不同。
不论是动态解还是静态解，在模型、源点深度和场点深度确定的情况下，格林函数仅与震中距相关。
实际建立静态格林函数库时，通常使用一维震中距序列；Python 接口中对应
``dists``，C 模块中对应 **-R**。结果在 NC 文件中保存为 ``north=0``、
``east=dists`` 的二维坐标，便于后续合成时在各震中距采样点上分别完成合成，
再按目标震中距加权组合合成结果。

模块也支持使用 ``norths`` / ``easts`` （C 模块 **-X/-Y**）指定二维北向/东向网格。
这种方式主要用于准确度测试，或希望后续合成直接复用同一网格、避免按震中距对多个合成结果进行加权组合的场景。
每个节点的震中距为 :math:`r_{ij} = \sqrt{north_i^2 + east_j^2}`，
实际计算中会对相同震中距自动去重以减少计算量。

结果会写入用户在构造 :class:`~pygrt.pymod.PyModel1D` 时通过 ``stgrn=`` 指定的 NetCDF 文件。
在 :doc:`static_syn` 阶段，可以指定新的网格或任意接收点，
程序会先利用格林函数库中包围目标震中距的采样点分别完成合成，再按目标震中距对这些合成结果进行加权组合；
需要根据深度进行权重组合时，对周围深度采样点的合成结果进行同样的加权组合。
使用与库一致的二维网格时可直接使用对应采样点的合成结果。

示例程序
-----------

假设在 :file:`milrow` 模型中，震源深度2km，接收点位于地表。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

        结果输出为 |NetCDF| 网格格式，方便使用 |GMT| 等软件处理和绘制。
        如下使用 ``ncdump -h`` 命令可查看网格文件基本信息。

        .. literalinclude:: run/grn_head
            :language: text


    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

        结果写入构造 :class:`~pygrt.pymod.PyModel1D` 时 ``stgrn=`` 指定的 NetCDF 文件。需要读回时调用
        :func:`pygrt.utils.read_static_nc`，返回字典包含
        ``dimensions``、 ``variables`` 与 ``attributes``。
