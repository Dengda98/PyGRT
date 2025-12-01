:author: 朱邓达
:date: 2025-04-24
:updated_date: 2025-09-22

固定间隔的Filon积分法
=========================

.. warning:: 

    谨慎使用固定间隔的Filon积分法，除非你很清楚原理以及误差来源。推荐使用 :doc:`/Advanced/filon/self_adaptive_filon` 以减少固定积分间隔的影响。 


以下介绍程序中使用的 **基于两点线性插值的Filon积分** :ref:`(纪晨, 姚振兴, 1995) <jichen_1995>`  :ref:`(初稿) <yao_init_manuscripts>` ，主要介绍思路，具体公式推导详见对应论文。

对于以下积分，

.. math:: 

   P_m(\omega,r) = \int_0^{k_{\text{max}}} F_m(k, \omega)J_m(kr)kdk 

分成两个部分，

.. math:: 

    P_m(\omega,r) = \int_0^{k^*} F_m(k, \omega)J_m(kr)kdk + \int_{k^*}^{k_{\text{max}}} F_m(k, \omega)J_m(kr)kdk 

其中积分第一部分仍然使用离散波数积分求解 （见 :doc:`/Advanced/k_integ` 部分），积分第二部分将Bessel函数取渐近表达式，可将其转为以下形式

.. math:: 

    \begin{align}
    P_m(\omega,r) &= \sqrt{\dfrac{2}{\pi r}}
    \bbox[yellow] {\int_{k^*}^{k_{\text{max}}} \bar{F}_m(k, \omega) \text{cos} \left( kr - (2m+1)\pi/4 \right) dk } \\
    \bar{F}_m(k, \omega) &= \sqrt{k} F_m(k, \omega)
    \end{align}

其中高亮部分的积分是Filon积分的关键。在每一个区间 :math:`[k_i, k_{i+1}], k_i=(i-1) \Delta k` 内，使用两点线性插值公式，:math:`\bar{F}` 可表示为 

.. math:: 

    \bar{F} = \bar{F}_i + \dfrac{\bar{F}_{i+1} - \bar{F}_i}{\Delta k} (k - k_i)

将此式代回积分公式中，可得到 **该区间内的积分解析解**，而最终的解即为多个区间内解析解的和。可进一步推导，将相邻区间的解析解中的部分式子两两相消，得到更简化的整个区间的解析解。由于相比于Bessel函数， :math:`\bar{F}` 往往变化比较缓慢，故Filon积分可使用相对较大的步长 :math:`\Delta k` 也可以有很好的精度。

------------------------------------

.. note:: 

    从以上推导可看出，该方法的精度取决于以下两个方面：

    1. **Bessel函数渐近表达式与真实值的误差：** 通过将积分分为两个部分计算，可大幅减少该误差。

    2. **区间内核函数线性拟合与真实值的误差：** 该误差受人为设定的Filon积分间隔控制。当 **区间过大** 时，线性拟合效果很差，计算结果偏差较大。


通过以下可选参数，可使用Filon积分。

.. tabs::

    .. group-tab:: C 

        详见 :doc:`/Module/greenfn` 和 :doc:`/Module/static_greenfn` 模块的 **-L** 选项。
         
    .. group-tab:: Python

        :func:`compute_grn() <pygrt.pymod.PyModel1D.compute_grn>` 函数和 :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` 函数支持以下可选参数来使用Filon积分，具体说明详见API。

        + ``Length:float``  定义离散波数积分的积分间隔 （见 :doc:`/Advanced/k_integ` 部分）
        + ``filonLength:float`` 定义Filon积分间隔，公式和离散波数积分使用的一致。
        + ``filonCut:float`` 定义了两个积分的分割点， :math:`k^*=` ``<filonCut>`` :math:`/r_{\text{max}}`

.. note::

    在程序中，固定间隔的Filon积分不计算近场项，即 :doc:`/Tutorial/dynamic/gfunc` 中介绍的 :math:`p=1` 对应的积分项。