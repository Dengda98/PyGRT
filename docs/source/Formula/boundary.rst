:author: 朱邓达
:date: 2026-01-07

自由边界和刚性边界的反射系数
============================================

本质上与  |yao2026| 中推导顶层的自由界面反射系数的过程基本一致，只是替换了相应的矩阵表达式。

具体而言， 以动态 P-SV 为例，

+ 自由边界满足过 z 平面的牵引力为 0， 即

.. math::

    \begin{bmatrix}
    q_m \\
    w_m \\
    0 \\
    0 \\
    \end{bmatrix} = 
    \begin{bmatrix}
    k   &   b   &   k  &  -b \\
    a   &   k   &  -a  &  k  \\
    2\mu\Omega & 2k\mu b & 2\mu\Omega & -2k\mu b \\
    2k\mu a & 2\mu\Omega & -2k\mu a & 2\mu\Omega \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^- \\
    \psi_m^- \\
    \phi_m^+ \\
    \psi_m^+ \\
    \end{bmatrix}

+ 刚性边界满足位移为 0 ，即

.. math::

    \begin{bmatrix}
    0 \\
    0 \\
    \sigma_{Rm} \\
    \tau_{Rm} \\
    \end{bmatrix} = 
    \begin{bmatrix}
    k   &   b   &   k  &  -b \\
    a   &   k   &  -a  &  k  \\
    2\mu\Omega & 2k\mu b & 2\mu\Omega & -2k\mu b \\
    2k\mu a & 2\mu\Omega & -2k\mu a & 2\mu\Omega \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^- \\
    \psi_m^- \\
    \phi_m^+ \\
    \psi_m^+ \\
    \end{bmatrix}

根据以上结果为 0 的表达式，可以构建出以下关系：

+ 对于顶层界面的 z+ 侧

.. math::

    \begin{bmatrix}
    \phi_m^+ \\
    \psi_m^+ \\
    \end{bmatrix} = 
    \mathbf{R}_U 
    \begin{bmatrix}
    \phi_m^- \\
    \psi_m^- \\
    \end{bmatrix}

+ 对于底层界面的 z- 侧

.. math::

    \begin{bmatrix}
    \phi_m^- \\
    \psi_m^- \\
    \end{bmatrix} = 
    \mathbf{R}_D 
    \begin{bmatrix}
    \phi_m^+ \\
    \psi_m^+ \\
    \end{bmatrix}


具体公式我使用 `SymPy <https://www.sympy.org/>`_ 做了推导，相应的 ``.ipynb`` 文件如下，以供读者参考。

:download:`boundary_condition.ipynb` （ :doc:`预览 <boundary_condition>` ）
