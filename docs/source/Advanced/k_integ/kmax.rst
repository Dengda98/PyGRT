:author: 朱邓达
:date: 2025-04-24

控制波数积分上限
======================

在 :doc:`/Tutorial/dynamic/gfunc` 部分介绍了程序将按以下形式从数值上计算波数积分，

.. math:: 

   P_m(\omega) = \Delta k \sum_{j=0}^{\infty} F_m(k_j,\omega)J_m(k_j r)k_j

其中 :math:`\Delta k = 2\pi/L, k_j=j\Delta k`，:math:`L` 为特征长度，即 :math:`L` 控制波数积分的积分间隔，默认根据 :doc:`/Tutorial/dynamic/gfunc` 部分介绍的约束条件自动确定。

程序会取一个较大值 :math:`k_{\text{max}}` 作为波数积分的上限，经验公式为

.. math:: 

    k_{\text{max}} = \sqrt{k_0^2 + \left(s \cdot \dfrac{\omega}{v_{\text{min}}}\right)^2}

其中 

+ :math:`k_0` 为零频的波数积分上限，默认为 :math:`\dfrac{5\pi}{h_s}`， :math:`h_s` 为震源和场点的深度差（绝对值），限制最小为1km，若实际深度差小于该值，则自动使用直接收敛法（DCM）；
+ :math:`\omega` 为角频率；
+ :math:`v_{\text{min}}` 为参考最小速度，默认取自模型中的最小速度，且限制在0.1km/s以上。
+ :math:`s` 为放大系数，默认为1.15；

程序支持自动判断积分收敛 |yao1983| 。当所有积分满足如下表达式时，自动退出波数循环（若达到 :math:`k_{\text{max}}` 则强制退出 ），

.. math:: 

    \left | \dfrac{ k_j F_m(k_j,\omega) J_m(k_j r) }
    {\sum_{i=1}^j k_i F_m(k_i,\omega) J_m(k_i r) } \right | \le \varepsilon

其中 :math:`\varepsilon` 为预先指定的收敛精度。程序默认不使用该功能。

------------------------------------

通过以下可选参数，可直接控制波数积分， **也将直接影响程序计算结果和运行时长**。

.. tabs:: 

    .. group-tab:: C 

        详见 :doc:`/Module/greenfn` 和 :doc:`/Module/static_greenfn` 模块的 **-K** 选项。

    .. group-tab:: Python

        :func:`compute_grn() <pygrt.pymod.PyModel1D.compute_grn>` 函数支持以下可选参数来控制波数积分，具体说明详见API。

        + ``k0:float``, 对应公式中的 :math:`k_0` 的系数，默认为5 
        + ``ampk:float``, 对应公式中的 :math:`s` ，默认为1.15 
        + ``keps:float`` 对应公式中的 :math:`\epsilon`，默认为-1（不使用）
        
        :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` 函数支持以下可选参数来控制波数积分，参数与上面对应，具体说明详见API。

        + ``k0:float`` 
        + ``keps:float`` 
