:author: 朱邓达
:date: 2025-04-17

计算静态应变、旋转、应力张量
=================================

除了使用不同的程序名和输出文件之外，流程基本和 :doc:`/Tutorial/dynamic/strain_stress` 类似。
Python 接口接收静态合成 NetCDF 路径，并由 C CLI 原地写入张量变量。
下面的示例先使用震中距序列建库，再在静态合成阶段指定二维接收网格，
结果展示水平面上的空间分布。


.. tabs:: 

    .. group-tab:: C 

        计算结果会以新增变量的形式直接写入 nc 网格，可使用 ``ncdump -h`` 查看。
        
        .. literalinclude:: run_upar/run.sh
            :language: bash
            :start-after: BEGIN
            :end-before: END


    .. group-tab:: Python

        .. literalinclude:: run_upar/run.py
            :language: python

-------------

.. figure:: run_upar/static_strain.svg
    :align: center 

-------------

.. figure:: run_upar/static_rotation.svg
    :align: center 

-------------

.. figure:: run_upar/static_stress.svg
    :align: center 


由于场点位于地表（自由表面），过Z平面的应力均为0（由于浮点数计算误差，呈极小非零数），结果和理论保持一致。
