:author: 朱邓达
:date: 2025-04-17
:updated_date: 2025-09-22

计算静态格林函数
=================

Python中计算静态格林函数的主函数为 :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` ，C模块为 :doc:`/Module/static_greenfn`。

建议先阅读完 :doc:`/Tutorial/dynamic/gfunc` 部分。静态情况与动态情况采取的计算方法一致，只是推导细节会有不同，详见 :ref:`初稿 <yao_init_manuscripts>` 。


示例程序
-----------

假设在 :file:`milrow` 模型中，震源深度2km，接收点位于地表。北向(X/km)在[-3,3]范围内等距间隔 0.15 km 采样，东向(Y/km)在[-2.5,2.5]范围内等距间隔 0.15 km 采样，计算这些点上的静态格林函数。

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