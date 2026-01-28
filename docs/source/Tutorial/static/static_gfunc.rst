:author: 朱邓达
:date: 2025-04-17

计算静态格林函数
=================

Python中计算静态格林函数的主函数为 :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` ，C模块为 :doc:`/Module/static_greenfn`。

建议先阅读完 :doc:`/Tutorial/dynamic/gfunc` 部分。静态情况与动态情况采取的计算方法一致，只是推导细节会有不同，详见 |yao2026| 。

静态解模块对于传入“震中距”的方式及后续处理与动态解模块相比稍有不同。
不论是动态解还是静态解，在模型、源点深度和场点深度确定的情况下，格林函数仅与 **和震中距** 相关。
但考虑到静态解的普遍应用场景，程序在计算静态格林函数阶段有两种传入震中距的方式：

+ 指定二维的 XY 网格。坐标 (0,0) 为源点的水平投影，每个节点的震中距为 :math:`r_{ij} = \sqrt{x_i^2 + y_j^2}` 。实际计算中会对震中距自动去重以减少计算量。
+ 指定一维的震中距序列。在计算和结果保存上等同于指定从原点 (X=0.0) 沿东向 (Y/km) 指定对应的采样点。

为了兼容性，结果的输出仍然保持以 XY 网格的形式。在 :doc:`static_syn` 阶段，可以指定新的 XY 网格，
此时每个节点的格林函数会近似为最近震中距的格林函数。

示例程序
-----------

假设在 :file:`milrow` 模型中，震源深度2km，接收点位于地表。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

        结果输出为 `NetCDF <https://zh.wikipedia.org/wiki/NetCDF>`_ 网格格式，方便使用 GMT 等软件处理和绘制。
        如下使用 ``ncdump -h`` 命令可查看网格文件基本信息。

        .. literalinclude:: run/grn_head
            :language: text


    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

        函数返回字典类型，包括一些基本参数以及格林函数（2D矩阵）。