:author: 朱邓达
:date: 2025-05-05

自适应Filon积分法
===================

在 :doc:`../linear_filon` 部分简单介绍了Filon积分的过程与核心原理，即 **分段计算积分的解析解**，而以下介绍的 **自适应Filon积分法（SAFIM）** ，便是 **动态分段**，在核函数变化剧烈的部分多采样，在核函数变化平缓的部分少采样。具体原理详见 :ref:`(Chen and Zhang, 2001) <chen_2001>` :ref:`(张海明, 2021) <zhang_book_2021>` 。


参数介绍
------------------

通过以下可选参数，可使用自适应Filon积分。

.. tabs:: 

    .. group-tab:: C 

        详见 :doc:`/Module/greenfn` 和 :doc:`/Module/static_greenfn` 模块的 **-L** 选项。

    .. group-tab:: Python

        :func:`compute_grn() <pygrt.pymod.PyModel1D.compute_grn>` 函数和 
        :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` 
        函数支持以下可选参数来使用自适应Filon积分，具体说明详见API。

        + ``Length:float``  定义离散波数积分的积分间隔 （见 :doc:`../kmax` 部分）
        + ``safilonTol:float`` 定义自适应采样精度，见 :ref:`(Chen and Zhang, 2001) <chen_2001>`  :ref:`(张海明, 2021) <zhang_book_2021>`，通常取 1e-2 即可。
        + ``filonCut:float`` 定义了两个积分的分割点， :math:`k^*=` ``<filonCut>`` :math:`/r_{\text{max}}` 


示例程序
------------------

.. tabs:: 

    .. group-tab:: C 

        使用 ``-S`` 导出核函数文件，

        .. literalinclude:: run/run.sh
            :language: bash 
            :start-after: BEGIN
            :end-before: END

    .. group-tab:: Python 

        .. literalinclude:: run/run.py
            :language: python 
    

使用以下Python脚本绘制核函数采样点，

.. literalinclude:: run/plot.py
    :language: python 


.. image:: run/safim.png
    :align: center


从对比图可以清晰地看出，自适应采样用更少的点数勾勒核函数的主要特征，这显著提高了计算效率。尽管自适应过程增加了额外的计算量，在小震中距范围内计算优势不明显， **但随着震中距的增加，自适应Filon积分的优势将愈发显著。**

.. note::

    在程序中，自适应Filon积分不计算近场项，即 :doc:`/Tutorial/dynamic/gfunc` 中介绍的 :math:`p=1` 对应的积分项。