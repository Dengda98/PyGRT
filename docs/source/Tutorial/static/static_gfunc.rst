:author: 朱邓达
:date: 2025-04-17

计算静态格林函数
=================

Python中计算静态格林函数的主函数为 :meth:`static_greenfn() <pygrt.pymod.PyModel1D.static_greenfn>` ，C模块为 :doc:`/Module/static_greenfn`。

建议先阅读完 :doc:`/Tutorial/dynamic/gfunc` 部分。静态情况与动态情况采取的计算方法一致，只是推导细节会有不同，详见 |yao2026p| 。


不论是动态解还是静态解，在模型、源点深度和场点深度确定的情况下，格林函数仅与震中距相关。
实际建立静态格林函数库时，通常直接指定一维震中距序列；Python 接口中对应
``dists``，C 模块中对应 **-R**。这是静态格林函数建库的常用方式。


示例程序
-----------

假设在 :file:`milrow` 模型中，震源深度2km，接收点位于地表。

.. tabs::  

    .. group-tab:: CLI

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

        结果输出为 |NetCDF| 文件。如下使用 ``ncdump -h`` 命令可查看文件的
        基本信息。输出文件仍保留 ``north`` 和 ``east`` 两个水平维度，这是为了兼容。

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
