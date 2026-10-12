:author: 朱邓达
:date: 2025-04-17

计算静态应变、旋转、应力张量
=================================

.. admonition:: 本页文档示例下载
    :class: tip

    :download:`static_strain_stress.tar.gz <static_strain_stress.tar.gz>` |
    解压后进入示例目录，运行 ``bash run.sh`` 即可复现本页示例。

除了使用不同的程序名和输出文件之外，流程基本和 :doc:`/Tutorial/dynamic/strain_stress/strain_stress` 类似。
计算结果直接写回输入的静态合成 NetCDF 文件。

**注意要在计算格林函数以及合成阶段的命令中加上计算位移偏导数的参数。**

以下示例中仅展示指定二维接收网格的用法，对任意接收点的用法也完全支持。


.. tabs:: 

    .. group-tab:: CLI

        计算结果会以新增变量的形式直接写入 nc 网格，可使用 ``ncdump -h`` 查看。
        
        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN
            :end-before: END


    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN
            :end-before: END

-------------

.. figure:: static_strain.svg
    :align: center 

-------------

.. figure:: static_rotation.svg
    :align: center 

-------------

.. figure:: static_stress.svg
    :align: center 


由于场点位于地表（自由表面），过Z平面的应力均为0（由于浮点数计算误差，呈极小非零数），结果和理论保持一致。


得到应力张量后，还可以进一步计算应力投影和库伦应力变化，
分别使用 :doc:`/Module/static_sproj` 和 :doc:`/Module/static_coulomb`。
投影由 **接收断层** 的走向、倾角和滑动角决定，与震源机制独立设置。
这里取接收断层走向 90°、倾角 60°、滑动角 0°，
将上面的应力张量投影为法向应力变化和沿滑动方向的剪应力变化。

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN SPROJ
            :end-before: END SPROJ

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN SPROJ
            :end-before: END SPROJ

结果 **sigma_n** 和 **tau_s** 写回 ``stsyn_dc_zne.nc``。
法向应力以张开为正，剪应力以沿接收断层滑动方向为正。
再取等效摩擦系数 0.4，计算库伦应力变化：

.. math::

    \Delta CFS = \Delta \tau_\mathrm{s} + \mu^{'} \Delta \sigma_\mathrm{n}

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN COULOMB
            :end-before: END COULOMB

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN COULOMB
            :end-before: END COULOMB

结果以 **coulomb** 变量写回同一文件，单位与应力相同，为 dyne/cm²（0.1 Pa）。


.. figure:: static_coulomb.svg
    :align: center

    从左到右为剪应力、法向应力和库伦应力变化，绘图单位为 MPa
