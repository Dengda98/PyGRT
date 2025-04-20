计算静态格林函数
=================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

Python中计算静态格林函数的主函数为 :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` ，C程序为 :command:`stgrt`。

建议先阅读完动态格林函数的计算过程（ :ref:`gfunc_rst` ）。静态情况与动态情况采取的计算方法一致，只是推导细节会有不同，详见 :ref:`初稿 <yao_init_manuscripts>` 。


示例程序
-----------

假设在 :file:`milrow` 模型中，震源深度2km，接收点位于地表。北向(X/km)在[-3,3]范围内等距采样41个点，东向(Y/km)在[-2.5,2.5]范围内等距采样33个点，计算这些点上的静态格林函数。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

        输出的 :file:`grn` 文件内部如下，开头的 ``#`` 用于保存源点和场点所在介质层的物性参数。

        .. literalinclude:: run/grn_head
            :language: text


    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

        函数返回字典类型，包括一些基本参数以及格林函数（2D矩阵）。