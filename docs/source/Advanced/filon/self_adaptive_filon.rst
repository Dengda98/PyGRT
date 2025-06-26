
.. _self_adaptive_filon:

自适应Filon积分法
===================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

在 :ref:`linear_filon` 部分简单介绍了Filon积分的过程与核心原理，即 **分段计算积分的解析解**，而以下介绍的 **自适应Filon积分法（SAFIM）** ，便是 **动态分段**，在核函数变化剧烈的部分多采样，在核函数变化平缓的部分少采样。具体原理详见 :ref:`(Chen and Zhang, 2001) <chen_2001>` :ref:`(张海明, 2021) <zhang_book_2021>` 。


参数介绍
------------------

通过以下可选参数，可使用自适应Filon积分。

.. tabs:: 

    .. group-tab:: C 

        :command:`grt` 和 :command:`stgrt` 程序支持以下可选参数来使用自适应Filon积分，具体说明详见 :command:`grt -h` 或 :command:`stgrt -h`。

        + ``-La<length>[/<Ftol>/<Fcut>]``
         
          + ``<length>``  定义离散波数积分的积分间隔 （见 :ref:`k_integ_rst` 部分, :ref:`linear_filon` 部分）
          + ``<Ftol>`` 定义自适应采样精度，见 :ref:`(Chen and Zhang, 2001) <chen_2001>`  :ref:`(张海明, 2021) <zhang_book_2021>`，通常1e-2即可。
          + ``<Fcut>`` 定义了两个积分的分割点， :math:`k^*=` ``<Fcut>`` :math:`/r_{\text{max}}` （见 :ref:`linear_filon` 部分）
         
    .. group-tab:: Python

        :func:`compute_grn() <pygrt.pymod.PyModel1D.compute_grn>` 函数和 :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` 函数支持以下可选参数来使用自适应Filon积分，具体说明详见API。

        + ``safilonTol:float`` 对应C选项卡中的 ``<Ftol>`` 参数
        + ``filonCut:float`` 对应C选项卡中的 ``<Fcut>`` 参数


示例程序
------------------

.. tabs:: 

    .. group-tab:: C 

        使用 ``-S`` 导出核函数文件，

        .. literalinclude:: run/run.sh
            :language: bash 

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