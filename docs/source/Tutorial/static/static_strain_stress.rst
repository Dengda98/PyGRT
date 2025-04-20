计算静态应变和应力
===================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

除了使用不同的函数名/程序名，输出文件不同之外，流程基本和 :ref:`strain_stress_rst` 类似。这里直接给出脚本。运行后注意观察和 :ref:`static_syn_rst` 的输出有了什么变化。

.. tabs:: 

    .. group-tab:: C 

        .. literalinclude:: run_upar/run.sh
            :language: bash


    .. group-tab:: Python

        .. literalinclude:: run_upar/run.py
            :language: python

-------------

.. image:: run_upar/static_strain.png
    :align: center 


-------------

.. image:: run_upar/static_stress.png
    :align: center 


由于场点位于地表（自由表面），过Z平面的应力均为0（由于浮点数计算误差，呈极小非零数），结果和理论保持一致。