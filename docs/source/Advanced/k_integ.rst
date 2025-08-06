控制波数积分
===================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

在 :doc:`/Tutorial/dynamic/gfunc` 部分介绍了程序将按以下形式从数值上计算波数积分，

.. math:: 

   P_m(\omega) = \Delta k \sum_{j=0}^{\infty} F_m(k_j,\omega)J_m(k_j r)k_j

其中 :math:`\Delta k = 2\pi/L, k_j=j\Delta k`，:math:`L` 为特征长度，即 :math:`L` 控制波数积分的积分间隔，默认根据 :doc:`/Tutorial/dynamic/gfunc` 部分介绍的约束条件自动确定。

程序会取一个较大值 :math:`k_{\text{max}}` 作为波数积分的上限，经验公式为

.. math:: 

    k_{\text{max}} = \sqrt{k_0^2 + s \cdot \left(\dfrac{\omega}{v_{\text{min}}}\right)^2}

其中 

+ :math:`k_0` 为零频的波数积分上限，默认为 :math:`\dfrac{5\pi}{h_s}`， :math:`h_s` 为震源和场点的深度差（绝对值），限制最小为1km，若实际深度差小于该值，则自动使用峰谷平均法；
+ :math:`\omega` 为角频率；
+ :math:`v_{\text{min}}` 为参考最小速度，默认取自模型中的最小速度，且限制在0.1km/s以上；若传入负值，表明使用峰谷平均法（负号不参与计算）。
+ :math:`s` 为放大系数，默认为1.15；

程序支持自动判断积分收敛 :ref:`(Yao and Harkrider, 1983) <yao&harkrider_1983>` 。当所有积分满足如下表达式时，自动退出波数循环（若达到 :math:`k_{\text{max}}` 则强制退出 ），

.. math:: 

    \left | \dfrac{ k_j F_m(k_j,\omega) J_m(k_j r) }
    {\sum_{i=1}^j k_i F_m(k_i,\omega) J_m(k_i r) } \right | \le \epsilon

其中 :math:`\epsilon` 为预先指定的收敛精度。程序默认不使用该功能。

------------------------------------

通过以下可选参数，可直接控制波数积分， **也将直接影响程序计算结果和运行时长**。

.. tabs:: 

    .. group-tab:: C 

        :command:`greenfn` 模块支持以下可选参数来控制波数积分，具体说明详见 :command:`grt greenfn -h`。  

        + ``-K<k0>[/<ampk>/<keps>]``  
         
          + ``<k0>`` 对应公式中的 :math:`k_0` 的系数，默认为5 
          + ``<ampk>`` 对应公式中的 :math:`s` ，默认为1.15
          + ``<keps>`` 对应公式中的 :math:`\epsilon`，默认为-1（不使用）
        
        + ``-V<vmin_ref>``， 对应公式中的 :math:`v_{\text{min}}`

        + ``-L<length>``， 对应公式中的 :math:`L` 的系数, :math:`L=` ``<length>`` :math:`\cdot r_{\text{max}}`
        
        :command:`static_greenfn` 模块支持以下可选参数来控制波数积分，参数与上面对应，具体说明详见 :command:`grt static greenfn -h`。

        + ``-K<k0>[/<keps>]``  
        + ``-L<length>``

    .. group-tab:: Python

        :func:`compute_grn() <pygrt.pymod.PyModel1D.compute_grn>` 函数支持以下可选参数来控制波数积分，具体说明详见API。

        + ``k0:float``, 对应公式中的 :math:`k_0` 的系数，默认为5 
        + ``ampk:float``, 对应公式中的 :math:`s` ，默认为1.15 
        + ``keps:float`` 对应公式中的 :math:`\epsilon`，默认为-1（不使用）
        + ``Length:float`` 对应公式中的 :math:`L` 的系数, :math:`L=` ``<length>`` :math:`\cdot r_{\text{max}}`
        
        :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` 函数支持以下可选参数来控制波数积分，参数与上面对应，具体说明详见API。

        + ``k0:float`` 
        + ``keps:float`` 
        + ``Length:float``  
