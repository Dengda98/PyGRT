:author: 朱邓达
:date: 2025-04-17

计算动态应变、旋转、应力张量
=======================================

计算这些张量需要计算位移的空间导数 :math:`\partial z,\partial r,\partial \theta` [#partial]_ 。 从前文关于格林函数的计算 （ :doc:`/Tutorial/dynamic/gfunc/gfunc` ） 和位移合成的公式 （ :doc:`/Tutorial/dynamic/syn/syn` ） 出发，位移空间导数分两个步骤计算，再根据几何方程和本构方程合成应变、旋转和应力张量。

.. [#partial] :math:`\partial z,\partial r,\partial \theta` 为 :math:`\dfrac{\partial}{\partial z},\dfrac{\partial}{\partial r},\dfrac{\partial}{\partial \theta}` 的简写。

计算格林函数阶段，求出 :math:`\partial z,\partial r`
-----------------------------------------------------------------------------------------------

假设在 :file:`milrow` 模型中，震源深度2km，接收点位于地表，震中距为30km，计算格林函数及其空间导数，要求采样点数1500个点，采样间隔0.02s。

.. tabs:: 

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

        在 :rst:dir:`GRN/milrow_2_0_30/` 路径下，文件名开头有"r"和"z"就代表 :math:`\partial r` 和 :math:`\partial z`。

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

        通道名开头有"r"和"z"就代表 :math:`\partial r` 和 :math:`\partial z`。


.. note:: 

    计算格林函数阶段求出的空间导数的单位：

    + 爆炸源：  :math:`10^{-25} \, /\text{dyne} \cdot \text{cm}`  
    + 单力源：  :math:`10^{-20} \, /\text{dyne}`
    + 剪切源：  :math:`10^{-25} \, /\text{dyne} \cdot \text{cm}`  
    + 矩张量源： :math:`10^{-25} \, /\text{dyne} \cdot \text{cm}`  



合成位移阶段，求出 :math:`\partial \theta / r` 并合成位移空间导数
-----------------------------------------------------------------------------

假设震源为剪切源，断层走向33°，倾角50°，滑动角120°，标量矩 1e24 dyne·cm，方位角30°。
本节使用阶跃源，并卷积上升 0.5 s、下降 1.5 s 的非对称余弦波。

.. tabs:: 

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN SYN DC
            :end-before: END SYN DC

        在 :rst:dir:`syn_dc/` 路径下，文件名开头有"r"，"z"，"t"分别代表 :math:`\partial r`， :math:`\partial z`，:math:`\partial \theta / r`。

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN SYN DC
            :end-before: END SYN DC

        通道名开头有"r"，"z"，"t"分别代表 :math:`\partial r`， :math:`\partial z`，:math:`\partial \theta / r`。


以上计算的空间导数在柱坐标系下，角向项使用 :math:`r^{-1}\partial_\theta`。
如果需要直角坐标系中的空间导数，可选择 ZNE 输出
（:math:`x` 表示北向，:math:`y` 表示东向），
此时文件名和通道名的前缀从 **z/r/t** 变为 **z/n/e**。

.. tabs:: 

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN ZNE
            :end-before: END ZNE

        在 :rst:dir:`syn_dc_zne/` 路径下，文件名开头有"z"，"n"，"e"分别代表 :math:`\partial z`， :math:`\partial x`，:math:`\partial y`。

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN ZNE
            :end-before: END ZNE

        通道名开头有"z"，"n"，"e"分别代表 :math:`\partial z`， :math:`\partial x`，:math:`\partial y`。

        .. warning::

            不可以使用ObsPy库中提供的Stream.rotate()函数进行位移空间导数的旋转。


.. note:: 

    合成的位移空间导数单位为1。


合成应变、旋转、应力张量
-----------------------------
应变和应力均为二阶对称张量，故将输出6个独立分量；而旋转张量为二阶反对称张量，将输出3个独立分量。

程序选择张量输出结果的坐标系的依据是——根据文件名/通道名判断是否存在 :math:`\dfrac{\partial u_x}{\partial x}`，如果有则使用ZNE坐标系，否则使用ZRT坐标系。 **所以建议保存结果的文件夹中只使用同种坐标系**，就像上面分为 :rst:dir:`syn_dc/` 和 :rst:dir:`syn_dc_zne/` 两个文件夹保存。

应变张量
~~~~~~~~~~~~~~
根据几何方程（应变-位移关系） [#covariant]_ 

.. math:: 

    e_{ij} = \dfrac{1}{2} \left( u_{i,j} + u_{j,i} \right) = \dfrac{1}{2} \left( \dfrac{\partial u_i}{\partial x_j} + \dfrac{\partial u_j}{\partial x_i}  \right)


.. tabs:: 

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN STRAIN
            :end-before: END STRAIN

        在 :rst:dir:`syn_dc_zne/` 原路径下，生成 :file:`strain_??.sac`，文件名中包括分量名，如ZZ、ZN等。

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN STRAIN
            :end-before: END STRAIN

        张量结果写回原合成目录中的 :file:`strain_??.sac`，
        通道名即为分量名，如 ZZ、ZN 等。
        动态结果只能靠文件名前缀区分应变/旋转/应力；
        示例脚本按 :file:`strain_*.sac` 读回，避免与位移等结果混在一起。

.. figure:: strain.svg
    :align: center


旋转张量
~~~~~~~~~~~~~~~~~~~~~
根据旋转-位移关系 [#covariant]_  

.. math:: 

    w_{ij} = \dfrac{1}{2} \left( u_{i,j} - u_{j,i} \right) = \dfrac{1}{2} \left( \dfrac{\partial u_i}{\partial x_j} - \dfrac{\partial u_j}{\partial x_i}  \right)

.. tabs:: 

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN ROTATION
            :end-before: END ROTATION

        在 :rst:dir:`syn_dc_zne/` 原路径下，生成 :file:`rotation_??.sac`，文件名中包括分量名，如ZN，ZE，NE。

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN ROTATION
            :end-before: END ROTATION

        张量结果写回原合成目录中的 :file:`rotation_??.sac`，
        通道名即为分量名，如 ZN、ZE、NE。
        示例脚本按 :file:`rotation_*.sac` 读回。

.. figure:: rotation.svg
    :align: center


应力张量
~~~~~~~~~~~~~~
根据各向同性介质的本构方程 [#covariant]_ 

.. math:: 

    \sigma_{ij} = \lambda \delta_{ij} e_{kk} + 2 \mu e_{ij} = \lambda \delta_{ij} u_{k,k} + \mu \left( u_{i,j} + u_{j,i} \right)



.. tabs:: 

    .. group-tab:: CLI

        .. literalinclude:: run.sh
            :language: bash
            :start-after: BEGIN STRESS
            :end-before: END STRESS

        在 :rst:dir:`syn_dc_zne/` 原路径下，生成 :file:`stress_??.sac`，文件名中包括分量名，如ZZ、ZN等。

    .. group-tab:: Python

        .. literalinclude:: run.py
            :language: python
            :start-after: BEGIN STRESS
            :end-before: END STRESS

        张量结果写回原合成目录中的 :file:`stress_??.sac`，
        通道名即为分量名，如 ZZ、ZN 等。
        示例脚本按 :file:`stress_*.sac` 读回。

.. figure:: stress.svg
    :align: center


.. note:: 

    计算得到的应力的单位为 :math:`\frac{\text{dyne}}{\text{cm}^2} = 0.1 \text{Pa}`。


场点位于地表（自由表面），**ZZ、ZN、ZE** 应力分量理论上为零。
图中这些分量使用整体应力的量级设置纵轴，避免自动缩放将微小数值残差放大。


得到应力张量后，还可以进一步计算应力投影和库伦应力变化，
分别使用 :doc:`/Module/sproj` 和 :doc:`/Module/coulomb`。
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

结果保存为 ``sigma_n.sac`` 和 ``tau_s.sac``。
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

结果保存为 ``coulomb.sac``，单位与应力相同，为 dyne/cm²（0.1 Pa）。

.. figure:: coulomb.svg
    :align: center

    从上到下为剪应力、法向应力和库伦应力变化，绘图单位为 MPa


.. [#covariant] 这只适用于ZNE坐标系，对于ZRT坐标系需考虑协变导数。程序中已考虑，这里不再做公式介绍。
