:author: 朱邓达
:date: 2026-06-24

直接收敛法的补充公式
============================================

在应用直接收敛法计算波数积分时，|dcm2025p| 附录给出了补偿项公式：

.. math::

    \int_0^\infty J_m(kr) k dk = 
    \left\{
    \begin{align}
    0, & \quad m = 0 \\
    \dfrac{1}{r^2}, & \quad m=1 \\
    \dfrac{2}{r^2}, & \quad m=2 \\
    \end{align}
    \right.
    

但以下情况需要拓展：

+ 计算位移的近场项（对应 :math:`p=1` 项的积分），详见 :doc:`/Tutorial/dynamic/gfunc` 
+ 计算位移对 :math:`z` 的偏导
+ 计算位移对 :math:`r` 的偏导

以下推导中需基于 |dcm2025p| 附录中提到的积分公式 (A3) 和 (A4)，以下将特定参数的积分结果展示如下，后续推导会使用到：

.. math:: 

    \int_0^\infty J_\nu(kr) k^\nu dk = 
    \left\{
    \begin{align}
    \dfrac{1}{r}, & \quad \nu = 0 \\
    \dfrac{1}{r^2}, & \quad \nu = 1 \\
    \dfrac{3}{r^3}, & \quad \nu = 2 \\
    \end{align}
    \right.

.. math::

    \int_0^\infty J_\nu(kr) k^{\nu+1} dk = 0, \quad \nu=0,1,2

**以下推导的公式均基于数值积分的方式进行了验证，详见** :download:`check_DCM_formula.ipynb` （ :doc:`预览 <check_DCM_formula>` ）

位移近场项
-----------------

位移的近场项涉及求解如下形式的积分，

.. math::

    \dfrac{1}{r} \int_0^\infty F(k) J_m(kr)dk

按照直接收敛法展开为

.. math::

    \dfrac{1}{r} \int_0^\infty F(k) J_m(kr)dk 
    = \dfrac{1}{r} \int_0^\infty \left[F(k) - F(k_\text{max}) \right] J_m(kr) dk + F(k_\text{max}) \cdot \bbox[yellow] {\dfrac{1}{r} \int_0^\infty J_m(kr) dk}

其中第二项即为直接收敛法的补偿项，其中高亮部分存在解析解，利用 Bessel 函数的性质：:math:`\int_0^\infty J_m(kr) dk = 1 / r` ，得到

.. math:: 

    \dfrac{1}{r} \int_0^\infty J_m(kr) dk = \dfrac{1}{r} \cdot \dfrac{1}{r} = \dfrac{1}{r^2}, \quad m=1,2

计算位移对 :math:`r` 的偏导
------------------------------------

位移对 :math:`r` 的偏导的效果将作用在 Bessel 函数上，以下将对远场项（对应 :math:`p=0,2,3` 项的积分）和近场项（对应 :math:`p=1` 项的积分）分别讨论。

-------------

对于 **远场项**，按照直接收敛法展开为

.. math::

    \int_0^\infty F(k) J_m^\prime(kr) k^2 dk = \int_0^\infty \left[F(k) - F(k_\text{max}) \right] J_m^\prime(kr) k^2 dk + F(k_\text{max}) \cdot \bbox[yellow] {\int_0^\infty J_m^\prime(kr) k^2 dk}

其中第二项的高亮部分存在解析解，对 :math:`m=0,1,2` 分别讨论：

+ :math:`m=0`

.. math::

    \int_0^\infty J_0^\prime(kr) k^2 dk = - \int_0^\infty J_1(kr) k^2 dk = - \lim\limits_{a \to 0^+} \dfrac{3ar}{(a^2 + r^2)^{5/2}} = 0

+ :math:`m=1`

.. math::

    & \quad \int_0^\infty J_1^\prime(kr) k^2 dk \\
    &= \int_0^\infty \left[ \dfrac{1}{kr} J_1(kr) - J_2(kr) \right] k^2 dk \\
    &= \dfrac{1}{r} \int_0^\infty J_1(kr) k dk - \int_0^\infty J_2(kr) k^2 dk \\
    &= \dfrac{1}{r} \cdot \dfrac{1}{r^2} - \dfrac{3}{r^3} \\
    &= - \dfrac{2}{r^3} 

+ :math:`m=2`

.. math::

    & \quad \int_0^\infty J_2^\prime(kr) k^2 dk \\
    &= \int_0^\infty \left[ J_1(kr) - \dfrac{2}{kr} J_2(kr) \right] k^2 dk \\
    &= \int_0^\infty J_1(kr) k^2 dk - \dfrac{2}{r} \int_0^\infty J_2(kr) k dk \\
    &= 0 - \dfrac{2}{r} \cdot \dfrac{2}{r^2} \\
    &=  - \dfrac{4}{r^3}

-------------

对于 **近场项**，按照直接收敛法展开为

.. math::
    \def \Br {\dfrac{\text{d}}{\text{d}r}\left[ \dfrac{1}{r} J_m(kr) \right]}

    & \quad \int_0^\infty F(k) \Br dk \\
    &= \int_0^\infty \left[F(k) - F(k_\text{max}) \right] \Br dk + F(k_\text{max}) \cdot \bbox[yellow] {\int_0^\infty \Br dk}

其中第二项的高亮部分存在解析解，进一步展开为

.. math::

    & \quad \int_0^\infty \dfrac{\text{d}}{\text{d}r}\left[ \dfrac{1}{r} J_m(kr) \right] dk \\
    &= \dfrac{1}{r} \int_0^\infty J_m^\prime(kr) k dk - \dfrac{1}{r^2} \int_0^\infty J_m(kr) dk \\
    &= \dfrac{1}{r} \int_0^\infty J_m^\prime(kr) k dk - \dfrac{1}{r^3}

对于上式的积分，对 :math:`m=1,2` 分别讨论：

+ :math:`m=1`

.. math::

    & \quad \int_0^\infty J_1^\prime(kr) k dk \\
    &= \int_0^\infty \left[ J_0(kr) - \dfrac{1}{kr} J_1(kr) \right] k dk \\
    &= \int_0^\infty J_0(kr) k dk - \dfrac{1}{r} \int_0^\infty J_1(kr) dk \\
    &= 0 - \dfrac{1}{r} \cdot \dfrac{1}{r}  \\
    &= - \dfrac{1}{r^2}

因此，

.. math:: 

    \dfrac{1}{r} \int_0^\infty J_1^\prime(kr) k dk - \dfrac{1}{r^3} = - \dfrac{2}{r^3} 

+ :math:`m=2`

.. math::

    & \quad \int_0^\infty J_2^\prime(kr) k dk \\
    &= \int_0^\infty \left[ J_1(kr) - \dfrac{2}{kr} J_2(kr) \right] k dk \\
    &= \int_0^\infty J_1(kr) k dk - \dfrac{2}{r} \int_0^\infty J_2(kr) dk \\
    &= \dfrac{1}{r^2} - \dfrac{2}{r} \cdot \dfrac{1}{r}  \\
    &= - \dfrac{1}{r^2}

因此，

.. math:: 

    \dfrac{1}{r} \int_0^\infty J_1^\prime(kr) k dk - \dfrac{1}{r^3} = - \dfrac{2}{r^3} 

-----------------

综上所述，应用直接收敛法计算位移对 :math:`r` 的偏导时将涉及以下补偿项：

.. math::

    \int_0^\infty J_m^\prime(kr) k^2 dk = 
    \left\{
    \begin{align}
    0, & \quad m = 0 \\
    - \dfrac{2}{r^3}, & \quad m=1 \\
    - \dfrac{4}{r^3}, & \quad m=2 \\
    \end{align}
    \right.

.. math::

    \int_0^\infty \dfrac{\text{d}}{\text{d}r}\left[ \dfrac{1}{r} J_m(kr) \right] dk 
    = - \dfrac{2}{r^3}, \quad m=1,2


计算位移对 :math:`z` 的偏导
------------------------------------

位移对 :math:`z` 的偏导的效果将作用在核函数上，以下将对远场项（对应 :math:`p=0,2,3` 项的积分）和近场项（对应 :math:`p=1` 项的积分）分别讨论。

-------------

对于 **远场项**，根据核函数的衰减特征，波数积分更适合进行如下形式的分解：

.. math::

    \int_0^\infty F(k) J_m(kr) k dk 
    = \int_0^\infty k \left[\dfrac{F(k)}{k} - \dfrac{F(k_\text{max})}{k_\text{max}} \right] J_m(kr) k dk 
      + \dfrac{F(k_\text{max})}{k_\text{max}} \cdot \bbox[yellow] { \int_0^\infty J_m(kr) k^2 dk}

其中第二项的高亮部分存在解析解，对 :math:`m=0,1,2` 分别讨论：

+ :math:`m=0`

.. math::

    & \quad \int_0^\infty J_0(kr) k^2 dk\\
    &= \int_0^\infty \left[ \dfrac{2}{kr} J_1(kr) - J_2(kr) \right] k^2 dk \\
    &= \dfrac{2}{r} \int_0^\infty J_1(kr) k dk - \int_0^\infty J_2(kr) k^2 dk \\
    &= \dfrac{2}{r} \cdot \dfrac{1}{r^2} - \dfrac{3}{r^3} \\
    &= - \dfrac{1}{r^3}

+ :math:`m=1`

.. math::

    \int_0^\infty J_1(kr) k^2 dk = 0

+ :math:`m=2`

.. math::

    \int_0^\infty J_2(kr) k^2 dk = \dfrac{3}{r^3}

--------------

对于 **近场项**，和上述一致，将波数积分分解为

.. math::

    \dfrac{1}{r} \int_0^\infty F(k) J_m(kr)dk 
    = \dfrac{1}{r} \int_0^\infty k \left[\dfrac{F(k)}{k} - \dfrac{F(k_\text{max})}{k_\text{max}} \right] J_m(kr) dk 
    + \dfrac{F(k_\text{max})}{k_\text{max}} \cdot \bbox[yellow] { \dfrac{1}{r} \int_0^\infty J_m(kr) k dk}

其中第二项的高亮部分存在解析解，只需在 |dcm2025p| 附录的补偿项公式上乘 :math:`1/r` 即可。

-----------------

综上所述，应用直接收敛法计算位移对 :math:`z` 的偏导时将涉及以下补偿项：

.. math::

    \int_0^\infty J_m k^2 dk = 
    \left\{
    \begin{align}
    - \dfrac{1}{r^3}, & \quad m = 0 \\
    0, & \quad m=1 \\
    \dfrac{3}{r^3}, & \quad m=2 \\
    \end{align}
    \right.

.. math::

    \dfrac{1}{r} \int_0^\infty J_m(kr) k dk = 
    \left\{
    \begin{align}
    \dfrac{1}{r^3}, & \quad m=1 \\
    \dfrac{2}{r^3}, & \quad m=2 \\
    \end{align}
    \right.
