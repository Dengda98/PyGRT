求解位移对 :math:`z` 的偏导
==============================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

求解应力、应变、旋转张量时需要求解偏导 :math:`\dfrac{\partial \mathbf{u}}{\partial z}` ，其本质只需要修改最终表达式中的“接收函数矩阵”即可。
推导过程其实也是推导z平面牵引力表达式的一部分。为了讲清推导过程，这里把一些公式再写一遍。

在 :ref:`初稿 <yao_init_manuscripts>` 以及 :doc:`RT` 中已经提到，三分量位移可展开在柱面谐矢量坐标系下，

.. math::
    :label:

    \def \Bm{\mathbf{B}_m}
    \def \Cm{\mathbf{C}_m}
    \def \Pm{\mathbf{P}_m}

    \mathbf{u} &= \sum_{m=0}^{\infty} \int_0^{\infty} \left( q_m \Bm + i v_m \Cm + w_m \Pm  \right) k dk \\

根据势函数 :math:`(\phi, \psi, \chi)` 与位移的关系，以及柱面谐矢量的表达式，可以得到位移 :math:`\mathbf{u}` 各分量的表达式，

.. math::
    :label:

    q_m &= k \phi_m + \dfrac{d \psi_m}{d z} \\
    w_m &= \dfrac{d \phi_m}{d z} + k \psi_m \\
    v_m &= k \chi_m \\

上式两边对 :math:`z` 求偏导，并根据垂直波函数的性质，同样引入 :math:`\epsilon` 这个符号变量（:math:`z` 正则正，:math:`z` 负则负）

.. math::
    :label:

    \def \parz#1{\dfrac{\partial #1}{d z}}

    \parz {q_m} &= -\epsilon ak \phi_m + b^2 \psi_m \\
    \parz {w_m} &= a^2 \phi_m - \epsilon bk \psi_m \\
    \parz {v_m} &= -\epsilon bk \chi_m \\

从而得到位移对 :math:`z` 的偏导与垂直波函数之间的关系，

.. math::
    :label:

    \def \parz#1{\dfrac{\partial #1}{d z}}

    \begin{bmatrix}
    \parz {q_m}  \\
    \parz {w_m}
    \end{bmatrix} &= 
    \begin{bmatrix}
    ak  &   b^2  &  -ak   &  b^2 \\
    a^2  &  bk   &   a^2  &   -bk \\
    \end{bmatrix}
    \begin{bmatrix}
    \phi_m^- \\
    \psi_m^- \\
    \phi_m^+ \\
    \psi_m^+ \\
    \end{bmatrix} \\
    \parz {v_m} &= 
    \begin{bmatrix}
    bk  &   -bk  \\
    \end{bmatrix}
    \begin{bmatrix}
    \chi_m^- \\
    \chi_m^+ \\
    \end{bmatrix}

这样在计算位移对应的垂直波函数时，只需改用上式的接收函数矩阵，即可得到偏导 :math:`\dfrac{\partial \mathbf{u}}{\partial z}` 。