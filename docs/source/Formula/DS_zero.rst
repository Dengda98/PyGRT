证明：当源点和场点均位于地表时，DS分量恒为0
============================================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

动态解
-------------

为了方便推导，我们假设 **场点位于源点下方** ，对应的 P-SV 位移核函数表达式为

.. math::
    :label: full_qw

    \def \Rm{\mathbf{R}}
    \def \Im{\mathbf{I}}
    \def \Tm{\mathbf{T}}

    \begin{bmatrix}
    q_m \\
    w_m \\
    \end{bmatrix}
    = 
    \Rm_{EV} \left( \Im - \Rm_U^{SR} \Rm_D^{RL} \right)^{-1} \Tm_D^{SR} \left( \Im - \Rm_U^{FS} \Rm_D^{SL} \right)^{-1}
    \bbox[yellow] {
    \begin{bmatrix}
        \begin{pmatrix}
        P_m^+ \\
        SV_m^+ \\
        \end{pmatrix}
        +
        \Rm_U^{FS}
        \begin{pmatrix}
        P_m^- \\
        SV_m^- \\
        \end{pmatrix}
    \end{bmatrix}
    }
    
其中高亮部分是关键。由于场点和源点都位于地表，矩阵 :math:`\mathbf{R}_U^{FS}` 退化为自由表面的反射系数矩阵 :math:`\tilde{\mathbf{R}}` ，即

.. math::
    :label: free_R

    \def \Rm{\mathbf{R}}

    \Rm_U^{FS} \rightarrow \tilde{\Rm} = \dfrac{1}{\Delta}
    \begin{bmatrix}
    k^2ab + \Omega^2  &  2kb\Omega \\
    2ka\Omega         &  k^2ab + \Omega^2 \\
    \end{bmatrix}

其中 :math:`\Delta=k^2ab-\Omega^2` ，:math:`\Omega = k^2 - k_\beta^2 / 2` 。

对于 P-SV 波， DS 分量的震源系数为

.. math::
    :label: coef_PSV

    P_1 &= 2 \varepsilon k \\
    SV_1 &= \dfrac{2k^2 - k_\beta^2}{b} \\

其中 :math:`\varepsilon` 为符号变量（对应 :eq:`full_qw` 式中震源系数的上标正负号）。
只需将 :eq:`coef_PSV` 式和 :eq:`free_R` 式代入 :eq:`full_qw` 式，即可证明 :eq:`full_qw` 式中高亮部分恒为0。

SH 波也有相同结论，此时自由界面的反射系数 :math:`\tilde{R}_L = 1` ，DS 分量的震源系数为 :math:`\chi_1 = \dfrac{-\varepsilon k_\beta^2}{k}` 。


静态解
-------------

证明过程和结论与动态解完全相同，这里列出证明会用到的自由界面反射系数矩阵和震源系数，读者可轻松证明。

自由界面的反射系数矩阵，其中 :math:`\Delta = \dfrac{\lambda + \mu}{\lambda + 3\mu}` ：

.. math::
    :label:

    \def \Rm{\mathbf{R}}

    \tilde{\Rm} &= 
    \begin{bmatrix}
    0                    &    - \Delta \\
    - \dfrac{1}{\Delta}  &    0  \\
    \end{bmatrix} \\
    \tilde{R}_L &= 1

DS 分量的震源系数：

.. math::
    :label:

    \def \eps{\varepsilon}
    \def \D{\Delta}

    P_1 &= - \dfrac{\eps \D}{1 + \D} \\
    SV_1 &= \dfrac{\eps }{1 + \D} \\
    SH_1 &= \eps \\
