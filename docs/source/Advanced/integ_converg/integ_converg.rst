:author: 朱邓达
:date: 2025-04-19

积分收敛性
===============================

通过输出核函数文件，观察源点和场点深度接近时，积分收敛性的变化。
**具体积分表达式以及分类详见** :doc:`/Tutorial/dynamic/gfunc`。

核函数文件
---------------
该文件记录了波数积分过程中不同波数对应的核函数值。以下用法展示了如何输出该文件。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN DGRN
            :end-before: END DGRN

        输出的核函数文件会在  :rst:dir:`GRN_grtstats/milrow_{depsrc}_{deprcv}/` 路径下。

    .. group-tab:: Python

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN DGRN
            :end-before: END DGRN

        输出的核函数文件会在自定义路径下。


C和Python导出的核函数文件是一致的，底层调用的是相同的函数。文件名称格式为 ``K_{iw}_{freq}``，其中 ``{iw}`` 表示频率索引值， ``{freq}`` 表示对应频率(Hz)。文件为自定义的二进制文件， **强烈建议使用Python进行读取及后续处理**。这里还是给出两种读取方法。

.. tabs::  

    .. group-tab:: C 

        :doc:`/Module/ker2asc` 模块可将单个核函数文件转为文本格式。

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN grt.ker2asc
            :end-before: END grt.ker2asc

        输出的文件如下，

        .. literalinclude:: run/stats_head
            :language: text

        后续你可以选择习惯的方式读取和处理。

    .. group-tab:: Python

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN read statsfile
            :end-before: END read statsfile

其中除了波数 ``k`` 外，每条结果的命名格式均为 ``{srcType}_{q/w/v}``，与 :doc:`/Tutorial/dynamic/gfunc` 部分介绍的积分公式中的核函数 :math:`q_m, w_m, v_m` 保持一致。

.. note:: 

    核函数文件中记录的值非理论核函数值。对于动态解，还需乘 :math:`\left(-\dfrac{1}{4\pi\rho\omega^2}\right)`。

可视化
-------------
以下将使用Python进行图件绘制。 **在Python函数中指定震源类型、阶数、积分类型，可自动绘制核函数、被积函数和积分值随波数的变化**，其中积分类型对应 :doc:`/Tutorial/dynamic/gfunc` 部分介绍的4种类型。

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN plot stats
    :end-before: END plot stats

.. image:: run/SS_0.png
    :align: center


可以通过指定 ``RorI=False`` 参数来指定绘制虚部，传入 ``RorI=2`` 表示实虚部都绘制。

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN plot stats RI
    :end-before: END plot stats RI

.. image:: run/SS_0_RI.png
    :align: center


从图中可以看到，当波数增加到一定范围以上，积分收敛良好，无振荡现象。



当震源和场点深度接近/相等时，积分收敛速度变慢
----------------------------------------------

假设其它参数不变，我们调整 **震源深度为0.1km**，再计算格林函数，读取核函数文件，绘制图像。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN DEPSRC 0.0 DGRN
            :end-before: END DEPSRC 0.0 DGRN

    .. group-tab:: Python

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN DEPSRC 0.0 DGRN
            :end-before: END DEPSRC 0.0 DGRN


.. image:: run/SS_0_0.0_RI.png
    :align: center

从图中可以清晰地看到，相比震源深度2km时，积分收敛速度明显变慢，积分值振荡，这要求增加波数积分上限，但这必然降低计算效率。

你可以尝试 **当震源和场点位于同一深度**，收敛振荡现象也很明显。要注意的是这和是否在地表无关，即使震源和场点都在地下，深度接近时积分收敛也会变慢。
