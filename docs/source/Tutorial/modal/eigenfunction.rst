:author: 朱邓达
:date: 2026-05-20

计算本征函数、群速度、频散敏感核
============================================

本征值对应的位移应力解被称为本征函数，再基于变分理论计算出能量积分后可计算出群速度和频散敏感核 (|aki2009|; |yao2026|)。

C 模块为 :doc:`/Module/eigenfn`，Python 中对应 :meth:`PyModel1D.eigenfn() <pygrt.pymod.PyModel1D.eigenfn>`，
计算将基于 :doc:`/Module/eigenv` 模块的频散结果及其内嵌的模型。

以下示例基于 :doc:`eigenvalue` 部分计算的相速度频散结果。

本征函数
------------------

计算频率 1 Hz 时不同阶的本征函数，深度范围 [0, 50] km，间隔 0.1 km :

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN EIGENFN
            :end-before: END EIGENFN

    .. group-tab:: Python

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN EIGENFN
            :end-before: END EIGENFN

        Python 提供了 :func:`read_eigenfunction() <pygrt.utils.read_eigenfunction>`
        将结果读入字典（``eigfn`` 等按阶号分组），以及 :func:`plot_eigenfunction() <pygrt.utils.plot_eigenfunction>`
        进行简易绘图。


:download:`plot_eigenfunction.py <run/plot_eigenfunction.py>`

.. figure:: run/eigenfunction.svg

群速度
-----------------

计算所有阶的相速度对应的群速度 :

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN GROUP
            :end-before: END GROUP

        :doc:`/Module/disp2asc` 模块也可将 |NetCDF| 格式的群速度频散转为文本格式。

    .. group-tab:: Python

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN GROUP
            :end-before: END GROUP

        与 :doc:`eigenvalue` 中介绍的一致，可使用
        :func:`read_dispersion() <pygrt.utils.read_dispersion>` 和
        :func:`plot_dispersion() <pygrt.utils.plot_dispersion>` 读取群速度并绘制；

以下同样提供两种方式的简易绘图脚本。

:download:`plot_dispersion_group.py (基于 Matplotlib) <run/plot_dispersion_group.py>` | 
:download:`plot_dispersion_group.sh (基于 GMT) <run/plot_dispersion_group.sh>`

.. figure:: run/dispersion_py_group.svg


频散敏感核
-----------------

计算不同频率下基阶频散曲线的相/群速度频散敏感核，深度间隔 0.2 km :

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SENS
            :end-before: END SENS

    .. group-tab:: Python

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SENS
            :end-before: END SENS

        Python 提供了 :func:`read_sensitivity() <pygrt.utils.read_sensitivity>`
        将结果读入字典（``csens`` / ``usens`` 等按阶号分组），以及 :func:`plot_sensitivity() <pygrt.utils.plot_sensitivity>`
        进行简易绘图。

:download:`plot_sensitivity.py <run/plot_sensitivity.py>`

+ **相速度频散敏感核**

.. grid:: 2

    .. grid-item::

        .. figure:: run/csens_R.svg
    
    .. grid-item::

        .. figure:: run/csens_L.svg


+ **群速度频散敏感核**

.. grid:: 2

    .. grid-item::

        .. figure:: run/usens_R.svg
    
    .. grid-item::

        .. figure:: run/usens_L.svg

.. note:: 

    这里值得一提的是，根据群速度和相速度的关系，

    .. math::

        U = \dfrac{c^2}{c - \omega \dfrac{\text{d} c}{\text{d} \omega}}

    群速度敏感核表达式为（其中 :math:`X` 表示 :math:`\alpha/\beta/\rho` ）,

    .. math::

        \dfrac{\text{d} U}{\text{d} X} = \dfrac{U}{c} \left(2 - \dfrac{U}{c} \right) \dfrac{\text{d} c}{\text{d} X}
        + \omega \left( \dfrac{U}{c} \right)^2 \dfrac{\text{d}}{\text{d} \omega} \left(\dfrac{\text{d} c}{\text{d} X} \right)

    其中可以看到计算群速度敏感核需要相速度敏感核在频率上的差分，因此无法完全解析地得到群速度敏感核。
    若频率间隔不够小，群速度敏感核计算精度就会下降。
