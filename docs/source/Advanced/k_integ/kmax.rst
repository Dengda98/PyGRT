:author: 朱邓达
:date: 2025-04-24

控制波数积分上限
======================

在 :doc:`/Tutorial/dynamic/gfunc` 部分介绍了程序将按以下形式从数值上计算波数积分，

.. math:: 

   P_m(\omega) = \Delta k \sum_{j=0}^{\infty} F_m(k_j,\omega)J_m(k_j r)k_j

其中 :math:`\Delta k = 2\pi/L, k_j=j\Delta k`，:math:`L` 为特征长度，即 :math:`L` 控制波数积分的积分间隔，
默认根据 :doc:`/Tutorial/dynamic/gfunc` 部分介绍的约束条件自动确定。

程序首先根据经验公式确定波数积分搜索区间的上界 :math:`k_{\text{max,ref}}` 。对动态全波解，

.. math:: 

    k_{\text{max,ref}} = \sqrt{\left(k_0 \cdot \dfrac{\pi}{h_s}\right)^2 + \left(s \cdot \dfrac{\omega}{v_{\text{min}}}\right)^2}

其中 

+ :math:`k_0` 为零频项的系数，默认为 50，程序内部使用 :math:`k_0 \cdot \pi / h_s` ；
  :math:`h_s=\max(|z_s-z_r|, 0.1)` km 为震源和场点的深度差；
+ :math:`\omega` 为角频率；
+ :math:`v_{\text{min}}` 为参考最小速度，默认取自模型中的最小速度，
  且限制在 0.1 km/s 以上；
+ :math:`s` 为放大系数（``ampk``），默认为 2.0。

对静态解，:math:`k_{\text{max,ref}} = k_0 \cdot \pi / h_s` 。

默认情况下，程序在 :math:`[\Delta k, k_{\text{max,ref}}]` 内基于核函数振幅搜索
实际积分上限 :math:`k_{\text{max}}` 。
同深度时判断核函数是否逼近常数，异深度时判断振幅是否衰减至 0 。
若指定 ``use_kmax_ref=True`` （C 模块 **+f**），则直接使用 :math:`k_{\text{max,ref}}` 作为
:math:`k_{\text{max}}` 。

若振幅搜索达到 :math:`k_{\text{max,ref}}` 仍未满足收敛判据，
程序在默认（Auto）模式下将自动启用直接收敛法（DCM）处理积分收敛。
当震源与场点完全同深度时，也会自动使用 DCM 。

程序还支持提前判断积分收敛 (|yao1983|) 。
当所有积分满足如下表达式时，自动退出波数循环（若达到 :math:`k_{\text{max}}` 则强制退出 ）。

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

        :func:`compute_grn() <pygrt.pymod.PyModel1D.compute_grn>` 函数支持以下可选参数来控制波数积分，
        具体说明详见API。

        + ``k0:float``, 对应公式中零频项的系数 :math:`k_0` ，默认为 50 
        + ``ampk:float``, 对应公式中的 :math:`s` ，默认为 2.0 
        + ``keps:float`` 对应公式中的 :math:`\epsilon`，默认为 -1（不使用）
        + ``use_kmax_ref:bool`` 为 True 时直接使用 :math:`k_{\text{max,ref}}` 作为
          积分上限
        
        :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` 函数支持以下可选参数来控制波数积分，
        参数与上面对应，具体说明详见API。

        + ``k0:float`` 
        + ``keps:float`` 
        + ``use_kmax_ref:bool``
