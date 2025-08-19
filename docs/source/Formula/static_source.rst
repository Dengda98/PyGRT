静力学单力源的震源系数
===========================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

在 :ref:`《理论地震图及其应用（初稿）》 <yao_init_manuscripts>` 以及相关文献中，已给出静力学位错震源的震源系数，这里不再重复。以下补充对应单力源的震源系数。

公式回顾
-----------------
求解静力学震源系数的方法与动力学方法类似。在无外力情况下，静力学弹性波方程为

.. math:: 
    :label:

    (\lambda+2\mu) \nabla(\nabla \cdot \mathbf{u}) - \mu \nabla \times (\nabla \times \mathbf{u}) = 0

此时位移 :math:`\mathbf{u}` 可使用三个线性独立的基本解表示，

.. math::
    :label: u_pot3

    \mathbf{u} = \nabla \phi
     + \left\{ 
        2 \dfrac{\partial\psi}{\partial z} \mathbf{e}_z - 
        \left[
            1 + 2 \Delta (z-z_j)\dfrac{\partial}{\partial z}
        \right] \nabla \psi
    \right\}
     + \nabla \times (\chi \mathbf{e}_z)

其中 :math:`\Delta=\dfrac{\lambda+\mu}{\lambda+3\mu}` ，位移势函数 :math:`\phi, \psi, \chi` 满足Laplace方程，

.. math::
    :label:

    \nabla^2 \phi &= 0 \\
    \nabla^2 \psi &= 0 \\
    \nabla^2 \chi &= 0 \\

在柱坐标系下，三个势函数的解可表示为

.. math::
    :label: pot3

    \phi(r, \theta, z) &= \sum_{m=-\infty}^{\infty} e^{im\theta} \int_0^\infty \phi_m(k,z) J_m(kr) k dk 
    \\
    \psi(r, \theta, z) &= \sum_{m=-\infty}^{\infty} e^{im\theta} \int_0^\infty \psi_m(k,z) J_m(kr) k dk 
    \\
    \chi(r, \theta, z) &= \sum_{m=-\infty}^{\infty} e^{im\theta} \int_0^\infty \chi_m(k,z) J_m(kr) k dk 
    \\

其中垂直波函数

.. math::
    :label: pot_BDF

    \phi_m &= B_m e^{-k|z-z_j|} \\
    \psi_m &= D_m e^{-k|z-z_j|} \\
    \chi_m &= F_m e^{-k|z-z_j|} \\


为方便求解，使用柱面谐矢量(vector cylindrical harmonics) :math:`(\mathbf{B}_m, \mathbf{C}_m, \mathbf{P}_m)` 表示位移 :math:`\mathbf{u}` ，

.. math::
    :label: u_qwv

    \mathbf{u} = \sum_{m=-\infty}^{\infty} e^{im\theta} \int_0^\infty 
    ( q_m(k,z) \mathbf{B}_m + v_m(k,z) \mathbf{c}_m + w_m(k,z) \mathbf{P}_m ) k dk

其中，

.. math::
    :label:

    \mathbf{B}_m &= \left(\mathbf{e}_r \dfrac{\partial}{\partial kr} + \mathbf{e}_\theta \dfrac{1}{kr}\dfrac{\partial}{\partial\theta}\right) J_m(kr)e^{im\theta} 
    \\
    \mathbf{C}_m &= \left(\mathbf{e}_r \dfrac{1}{kr}\dfrac{\partial}{\partial\theta} - \mathbf{e}_\theta \dfrac{\partial}{\partial kr}\right) J_m(kr)e^{im\theta} 
    \\
    \mathbf{P}_m &= \mathbf{e}_z J_m(kr)e^{im\theta}

将 :eq:`u_pot3` 式和 :eq:`pot3` 式代入 :eq:`u_qwv` 式，得到系数 :math:`q_m, v_m, w_m` 和垂直波函数之间的关系，

.. math:: 
    :label: qwv_pot

    q_m &= \phi_m - \left[ 1 + 2\Delta(z-z_j)\dfrac{\partial}{\partial z}  \right] \psi_m 
    \\
    w_m &= \dfrac{1}{k} \dfrac{\partial\phi_m}{\partial z} + \dfrac{1}{k} \left[ \dfrac{\partial}{\partial z} - 2\Delta k^2 (z-z_j)  \right] \psi_m 
    \\
    v_m &= \chi_m 
    \\

其中三个垂直波函数均吸收了 :math:`k` 因子，即 :math:`\phi_m \leftarrow k\phi_m, \psi_m \leftarrow k\psi_m, \chi_m \leftarrow k\chi_m` ，这会体现在后续的震源系数中。

为了求解震源系数，需将 :doc:`static_uniform` 部分的最终表达式展开成 :eq:`u_qwv` 的形式，对应的垂直波函数的系数 :math:`(B_m, D_m, F_m)` 即为震源系数。具体而言，我们可将位移同一分量的不同表达式进行对比得到震源系数。将 :eq:`qwv_pot` 式代入 :eq:`u_qwv` 式，柱坐标系下的位移三分量表达式为

.. math::
    :label:

    u_r &= \sum_{m=-\infty}^{\infty} e^{im\theta} \int_0^\infty 
    \left\{
    \left\{
        \phi_m -  \left[ 1 + 2\Delta(z-z_j) \dfrac{\partial}{\partial z} \right] \psi_m
    \right\} J_m^{'}(kr)
    +  \chi_m \dfrac{im}{kr} J_m(kr)
    \right\} k dk 
    \\
    u_\theta &= \sum_{m=-\infty}^{\infty} i e^{im\theta} \int_0^\infty 
    \left\{
    \left\{
         \phi_m -  \left[ 1 + 2\Delta(z-z_j) \dfrac{\partial}{\partial z} \right] \psi_m
    \right\} \dfrac{m}{kr} J_m(kr)
    + i \chi_m J_m^{'}(kr)
    \right\} k dk 
    \\
    u_z &= \sum_{m=-\infty}^{\infty} e^{im\theta} \int_0^\infty 
    \left\{
        \dfrac{1}{k} \dfrac{\partial\phi_m}{\partial z} + \dfrac{1}{k} \left[ \dfrac{\partial}{\partial z} - 2\Delta k^2 (z-z_j) \right] \psi_m
    \right\} J_m(kr) k dk 
    \\
    
再将 :eq:`pot_BDF` 式代入，整理得到，

.. math::
    :label: u_BDF

    u_r &= \sum_{m=-\infty}^{\infty} e^{im\theta} \int_0^\infty 
    \left\{
    \left[
        (B_m - D_m) + 2\Delta\varepsilon k  (z-z_j) D_m
    \right] J_m^{'}(kr)
    +  F_m \dfrac{im}{kr} J_m(kr)
    \right\} e^{-k|z-z_j|} k dk 
    \\
    u_\theta &= \sum_{m=-\infty}^{\infty} i e^{im\theta} \int_0^\infty 
    \left\{
    \left[
        (B_m - D_m) + 2\Delta\varepsilon k  (z-z_j) D_m
    \right] \dfrac{m}{kr} J_m(kr)
    + i  F_m J_m^{'}(kr)
    \right\} e^{-k|z-z_j|} k dk 
    \\
    u_z &= \sum_{m=-\infty}^{\infty} e^{im\theta} \int_0^\infty 
    \left[-\varepsilon  (B_m + D_m) - 2\Delta k (z-z_j) D_m \right]
    J_m(kr) e^{-k|z-z_j|} k dk

其中，

.. math::

    \varepsilon = 
    \begin{cases}
    -1, \ \ &当 z < z_j 为上行波 \\
    1,  \ \ & 当 z > z_j 为下行波 \\
    \end{cases}


垂向力源的震源系数
----------------------
当力的方向为 :math:`\mathbf{n}=\mathbf{e}_z` 时，径向的静态位移为

.. math::
    :label: vert_ur

    u_r = G_{13} \cos \theta + G_{23} \sin \theta = \dfrac{1}{4\pi \mu} \dfrac{\Delta}{1+\Delta} \dfrac{r(z-z_j)}{R^3} 

其中 :math:`r` 为水平距离， :math:`R=\sqrt{r^2 + (z-z_j)^2}` 为水平距离。根据Sommerfeld积分公式，有

.. math::
    :label:

    \dfrac{r(z-z_j)}{R^3} = \int_0^\infty  k (z-z_j) \dfrac{e^{-k|z-z_j|}}{k} J_1(kr) k dk 

将上式代入 :eq:`vert_ur` 式，利用Bessel函数的性质 :math:`J_0^{'}(x) = -J_1(x)` ，得到径向位移的积分表达式

.. math::
    :label:

    u_r = \dfrac{1}{4\pi \mu} \dfrac{\Delta}{1+\Delta} \int_0^\infty - (z-z_j) e^{-k|z-z_j|} J_0^{'}(kr) k dk 

将上式与 :eq:`u_BDF` 式中的径向位移 :math:`u_r` 表达式对比，得到震源系数

.. math::
    :label:

    B_0 &= D_0 = - \dfrac{\varepsilon}{8\pi\mu(1+\Delta)k} \\
    F_0 &= 0

为简化表达式和系数，将 :eq:`pot_BDF` 式代入 :eq:`u_pot3` 式，结合以上解的具体形式，提出公共因子，整理得到势函数的表达式，

.. math::
    :label: pot_A_PS

    \phi(r, \theta, z) &= \dfrac{1}{4\pi\mu} \sum_{m=0}^{2} A_m(\theta) \int_0^\infty P_m(k,z) e^{-k|z-z_j|}  J_m(kr) k dk 
    \\
    \psi(r, \theta, z) &= \dfrac{1}{4\pi\mu} \sum_{m=0}^{2} A_m(\theta) \int_0^\infty SV_m(k,z) e^{-k|z-z_j|}  J_m(kr) k dk 
    \\
    \chi(r, \theta, z) &= \dfrac{1}{4\pi\mu} \sum_{m=1}^{2} A_{m+2}(\theta) \int_0^\infty SH_m(k,z) e^{-k|z-z_j|}  J_m(kr) k dk 
    \\

其中方向因子为

.. math::
    :label:

    A_0 = 1

震源系数为

.. math::
    :label:

    P_0 &= SV_0 = - \dfrac{\varepsilon}{2(1+\Delta)k} \\
    SH_0 &= 0 


水平力源的震源系数
-------------------------

当力的方向为 :math:`\mathbf{n}=\mathbf{e}_x` 时，径向的静态位移为

.. math::
    :label: hori_ur

    u_r &= G_{11} \cos \theta + G_{21} \sin \theta \\
    &= \dfrac{1}{4\pi \mu} \dfrac{\Delta}{1+\Delta} \left( \dfrac{1}{\Delta R} \cos\theta + \dfrac{r^2\cos^2\theta}{R^3} \cos\theta + \dfrac{r^2\sin^2\theta}{R^3}  \cos\theta \right) \\
    &= \dfrac{1}{4\pi \mu} \dfrac{\Delta}{1+\Delta} \left( \dfrac{1}{\Delta R} + \dfrac{r^2}{R^3} \right) \cos\theta

根据Sommerfeld积分公式，有

.. math::
    :label: smfld_1R

    \dfrac{1}{R} &= \int_0^\infty \dfrac{e^{-k|z-z_j|}}{k}  J_0(kr) k dk \\
    \dfrac{r^2}{R^3} &= \int_0^\infty \left[ 1 - \varepsilon k (z-z_j) \right] \dfrac{e^{-k|z-z_j|}}{k} J_0(kr) k dk

利用Bessel函数的性质 :math:`J_0(x) = J_1^{'}(x) + \dfrac{1}{x} J_1(x)` ，得到

.. math::
    :label:

    \dfrac{1}{R} &= \int_0^\infty \left[ J_1^{'}(kr) + \dfrac{1}{kr} J_1(kr) \right] \dfrac{e^{-k|z-z_j|}}{k} k dk \\
    \dfrac{r^2}{R^3} &=  \int_0^\infty \left[ 1 - \varepsilon k (z-z_j) \right] \left[ J_1^{'}(kr) + \dfrac{1}{kr} J_1(kr) \right] \dfrac{e^{-k|z-z_j|}}{k} k dk \\
    &= \int_0^\infty 
        \left[ 
            -\varepsilon k (z-z_j) J_1^{'}(kr) + \underline{J_1^{'}(kr) -\varepsilon k (z-z_j) \dfrac{1}{kr} J_1(kr)} + \dfrac{1}{kr} J_1(kr) 
        \right]  \dfrac{e^{-k|z-z_j|}}{k} k dk 

注意下划线的两项，以下使用分部积分法证明其积分为0，

.. math::
    :label:

    &\int_0^\infty J_1^{'}(kr) e^{-k|z-z_j|} dk \\
    = &\int_0^\infty \dfrac{1}{r}  e^{-k|z-z_j|} d J_1(kr) \\
    = &\dfrac{1}{r} e^{-k|z-z_j|} J_1(kr) \Big|_0^\infty 
       + \int_0^\infty \varepsilon (z-z_j) \dfrac{1}{r} J_1(kr) e^{-k|z-z_j|} dk \\
    = &\int_0^\infty \varepsilon (z-z_j) \dfrac{1}{r} J_1(kr) e^{-k|z-z_j|} dk \\

因此，

.. math::
    :label: smfld_r2R3

    \dfrac{r^2}{R^3} = \int_0^\infty 
        \left[ 
            -\varepsilon k (z-z_j) J_1^{'}(kr) + \dfrac{1}{kr} J_1(kr) 
        \right]  \dfrac{e^{-k|z-z_j|}}{k} k dk 

将 :eq:`smfld_1R` 式和 :eq:`smfld_r2R3` 式代入 :eq:`hori_ur` 式，得到

.. math::
    :label:

    u_r = \dfrac{1}{4\pi \mu} \dfrac{\Delta}{1+\Delta} \cos\theta \int_0^\infty 
    \left\{
        \left[ \dfrac{1}{\Delta k} - \varepsilon (z-z_j) \right] J_1^{'}(kr) 
        + \dfrac{1+\Delta}{\Delta} \dfrac{1}{k} \dfrac{1}{kr} J_1(kr) 
    \right\} e^{-k|z-z_j|} k dk

将上式与 :eq:`u_BDF` 式中的径向位移 :math:`u_r` 表达式对比，得到震源系数

.. math::
    :label:

    B_1 &= - D_1 = \dfrac{1}{8\pi\mu(1+\Delta)k} \\
    F_1 &= - \dfrac{i}{4\pi \mu k}

同样将势函数表达成 :eq:`pot_A_PS` 式，得到方向因子 [#]_

.. math::
    :label:

    A_1 &= \cos\theta \\
    A_4 &= - \sin\theta

震源系数

.. math::
    :label:

    P_1 &= - SV_1 = \dfrac{1}{2 (1+\Delta)k} \\
    SH_1 &= - \dfrac{1}{k}


对比动态解的推导过程，静态解除了震源系数不同，方向因子均一致，后续的矩阵传播逻辑也一致，因此程序中动态解和静态解共享函数模块。

--------------------------------

.. [#] 水平力源的 :math:`A_4` 符号与 :ref:`《理论地震图及其应用（初稿）》 <yao_init_manuscripts>` 中所使用的相反，对应的方向因子 :math:`SH_m` 也相反，这对最终位移计算结果无影响。




