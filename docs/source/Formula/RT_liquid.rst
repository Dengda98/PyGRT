:author: 朱邓达
:date: 2025-06-23
:updated_date: 2025-08-19

含液体层的反射透射系数矩阵
==================================

当界面两侧存在液体时，边界条件变为 **垂直位移连续，垂直应力连续，切向应力为0** ，此时使用 :doc:`RT` 介绍的方法，可以比较方便直观地推导含液体层的 R/T 矩阵。

以下推导以动态解为例，最终每一项系数的详细表达式我使用 Python 库 `SymPy <https://www.sympy.org/>`_ 辅助推导，可在这里下载 :download:`RT_liquid_formula.ipynb` （ :doc:`预览 <RT_liquid_formula>` ）。

液体-液体界面
------------------

在液体层内，位移垂直波函数（仅P波）与位移应力矢量之间的关系可以表示为（略去角标 :math:`j`）

.. math::
    :label:

    \begin{bmatrix}
    w_m \\
    \sigma_{Rm} \\
    \end{bmatrix} = 
    \begin{bmatrix}
    a              &  -a             \\
    -\rho\omega^2  & -\rho\omega^2   \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^- \\
    \phi_m^+ \\
    \end{bmatrix}

相当于从固体情况的关系式中提取部分项，其中 :math:`2\mu\Omega` 转变为

.. math::
    
    2\mu\Omega \rightarrow 
    \lim_{\beta \rightarrow 0} 2\beta^2\rho(k^2 - \dfrac{1}{2}k_\beta^2)
    = -\rho\omega^2

为简单起见，我们假设界面位于第1/2层之间，根据边界条件，有

.. math::
    :label:

    \begin{bmatrix}
    a_1              &  -a_1             \\
    -\rho_1\omega^2  & -\rho_1\omega^2   \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^-(z_1^-) \\
    \phi_m^+(z_1^-) \\
    \end{bmatrix} = 
    \begin{bmatrix}
    a_2              &  -a_2             \\
    -\rho_2\omega^2  & -\rho_2\omega^2   \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^-(z_1^+) \\
    \phi_m^+(z_1^+) \\
    \end{bmatrix}

后续的 R/T 矩阵推导过程与固体-固体界面情况完全一样，只是矩阵元素和维度发生变化（但仍是方阵）。通过类似的对波入射方向的讨论，可以得到 R/T 矩阵满足（此时 R/T 矩阵已退化为 1x1 的标量）

.. math::
    :label:

    \begin{bmatrix}
    T_U & R_D \\
    R_U  & T_D \\
    \end{bmatrix} = 
    \begin{bmatrix}
    -a_1              &  -a_2             \\
    \rho_1\omega^2    & -\rho_2\omega^2   \\
    \end{bmatrix}^{-1}
    \begin{bmatrix}
    -a_2              &  -a_1             \\
    \rho_2\omega^2    & -\rho_1\omega^2   \\
    \end{bmatrix}

其中

.. math::
    :label:

    R_D &= \dfrac{a_1\rho_2 - a_2\rho_1}{a_1\rho_2 + a_2\rho_1} \\
    R_U &= \dfrac{a_2\rho_1 - a_1\rho_2}{a_1\rho_2 + a_2\rho_1} \\
    T_D &= \dfrac{2 a_1\rho_1}{a_1\rho_2 + a_2\rho_1} \\
    T_U &= \dfrac{2 a_2\rho_2}{a_1\rho_2 + a_2\rho_1} \\

在程序中为保持 2x2 矩阵，只将以上结果填充在对应位置即可，其余项为0，例如

.. math::
    
    \mathbf{R}_D^{2\times2} = 
    \begin{bmatrix}  
    R_D  &  0  \\
    0    &  0  \\
    \end{bmatrix}

液体-固体界面
------------------

我们使用角标 :math:`l` 表示液体层，角标 :math:`s` 表示固体层。当上层为液体，下层为固体时，根据边界条件，有

.. math::
    :label: layer12

    \left[
    \begin{array}{c|c}
    a_l              &  -a_l             \\
    -\rho_l\omega^2  & -\rho_l\omega^2   \\
    \hline
    0                &  0                \\
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \phi_m^-(z_1^-) \\
    \hline
    \phi_m^+(z_1^-) \\
    \end{array}
    \right] = 
    \left[
    \begin{array}{cc|cc}
    a_s   &   k     &  -a_s  &  k  \\
    2\mu_s\Omega_s & 2k\mu_s b_s & 2\mu_s\Omega_s & -2k\mu_s b_s \\
    \hline
    2k\mu_s a_s & 2\mu_s\Omega_s & -2k\mu_s a_s & 2\mu_s\Omega_s \\
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \phi_m^- (z_1^+) \\
    \psi_m^- (z_1^+) \\
    \hline
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{array}
    \right]

其中上层液体层仅有 P 波垂直波函数（2x1），而下层固体层有的 P、SV 波垂直波函数，这直接表明此情况下的 R/T 矩阵形状不一，会存在 2x1、1x2 等的形状。

波从上向下入射
~~~~~~~~~~~~~~~~

此时下层没有向上传播的入射波，即 :math:`[\phi_m^- (z_1^+), \psi_m^- (z_1^+)]^T = \mathbf{0}` ，:eq:`layer12` 式变为

.. math::
    :label:

    \left[
    \begin{array}{c|c}
    a_l              &  -a_l             \\
    -\rho_l\omega^2  & -\rho_l\omega^2   \\
    \hline
    0                &  0                \\
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \phi_m^-(z_1^-) \\
    \hline
    \bbox[yellow] {\phi_m^+(z_1^-)} \\
    \end{array}
    \right] = 
    \left[
    \begin{array}{cc}
    -a_s  &  k  \\
    2\mu_s\Omega_s & -2k\mu_s b_s \\
    \hline
    -2k\mu_s a_s & 2\mu_s\Omega_s \\
    \end{array}
    \right]
    \begin{bmatrix}
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{bmatrix}

其中高亮部分的垂直波函数为当前情况的“已知项”，通过移项+矩阵重排的方式可得到

.. math::
    :label: U2D

    \left[
    \begin{array}{c|cc}
    -a_l            &     -a_s           &  k             \\
    \rho_l\omega^2  &     2\mu_s\Omega_s & -2k\mu_s b_s   \\
    \hline
    0               &     -2k\mu_s a_s   & 2\mu_s\Omega_s \\
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \phi_m^-(z_1^-) \\
    \hline
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{array}
    \right] = 
    \left[
    \begin{array}{c}
    -a_l \\
    -\rho_l\omega^2 \\
    \hline
    0
    \end{array}
    \right]
    \left[\bbox[yellow]{\phi_m^+(z_1^-)}\right]

其中等号左边矩阵前两列的负号由移项产生，此时左边的垂直波函数矢量（作为未知量）已经变成两层的混合版本，适定方程可简单使用逆矩阵求解，得到

.. math::
    :label:

    \begin{bmatrix}
    \phi_m^- (z_1^-) \\
    \end{bmatrix} = 
    \mathbf{R}_D^{1\times1}
    \begin{bmatrix}
    \bbox[yellow] {\phi_m^+ (z_1^-)} \\
    \end{bmatrix} 

    \begin{bmatrix}
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{bmatrix} = 
    \mathbf{T}_D^{2\times1}
    \begin{bmatrix}
    \bbox[yellow] {\phi_m^+ (z_1^-)} \\
    \end{bmatrix} 


波从下向上入射
~~~~~~~~~~~~~~~~

此时上层没有向下传播的入射波，即 :math:`\phi_m^+ (z_1^-) = 0` ，:eq:`layer12` 式变为

.. math::
    :label:

    \left[
    \begin{array}{c}
    a_l \\
    -\rho_l\omega^2 \\
    \hline
    0
    \end{array}
    \right]
    \left[\phi_m^-(z_1^-)\right] = 
    \left[
    \begin{array}{cc|cc}
    a_s   &   k     &  -a_s  &  k  \\
    2\mu_s\Omega_s & 2k\mu_s b_s & 2\mu_s\Omega_s & -2k\mu_s b_s \\
    \hline
    2k\mu_s a_s & 2\mu_s\Omega_s & -2k\mu_s a_s & 2\mu_s\Omega_s \\
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \bbox[yellow] {\phi_m^- (z_1^+)} \\
    \bbox[yellow] {\psi_m^- (z_1^+)} \\
    \hline
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{array}
    \right]

其中高亮部分的垂直波函数同样为当前情况的“已知项”，为保持与 :eq:`U2D` 式的形式匹配，通过类似的移项+矩阵重排的方式可得到

.. math::
    :label: D2U

    \left[
    \begin{array}{c|cc}
    -a_l            &     -a_s           &  k             \\
    \rho_l\omega^2  &     2\mu_s\Omega_s & -2k\mu_s b_s   \\
    \hline
    0               &     -2k\mu_s a_s   & 2\mu_s\Omega_s \\
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \phi_m^-(z_1^-) \\
    \hline
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{array}
    \right] = 
    \left[
    \begin{array}{cc|cc}
    -a_s   &   -k    \\
    -2\mu_s\Omega_s & -2k\mu_s b_s \\
    \hline
    -2k\mu_s a_s & -2\mu_s\Omega_s  \\
    \end{array}
    \right]
    \left[
    \begin{array}{c}
    \bbox[yellow] {\phi_m^- (z_1^+)} \\
    \bbox[yellow] {\psi_m^- (z_1^+)} \\
    \end{array}
    \right]

矩阵中的负号由移项产生，等号左边形式与 :eq:`U2D` 式完全一致。同样该适定方程可简单使用逆矩阵求解，得到

.. math::
    :label:

    \begin{bmatrix}
    \phi_m^- (z_1^-) \\
    \end{bmatrix} = 
    \mathbf{T}_U^{1\times2}
    \begin{bmatrix}
    \bbox[yellow] {\phi_m^- (z_1^+)} \\
    \bbox[yellow] {\psi_m^- (z_1^+)} \\
    \end{bmatrix} 

    \begin{bmatrix}
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{bmatrix} = 
    \mathbf{R}_U^{2\times2}
    \begin{bmatrix}
    \bbox[yellow] {\phi_m^- (z_1^+)} \\
    \bbox[yellow] {\psi_m^- (z_1^+)} \\
    \end{bmatrix} 

合并求解
~~~~~~~~~~

:eq:`U2D` 式和 :eq:`D2U` 式可合并，一并使用逆矩阵求得最终液体-固体界面上的 R/T 矩阵，

.. math::
    :label:

    \begin{bmatrix}
    \mathbf{T}_U^{1\times2}  & \mathbf{R}_D^{1\times1} \\
    \mathbf{R}_U^{2\times2}  & \mathbf{T}_D^{2\times1} \\
    \end{bmatrix}_{3\times3} = 
    \left[
    \begin{array}{c|cc}
    -a_l            &     -a_s           &  k             \\
    \rho_l\omega^2  &     2\mu_s\Omega_s & -2k\mu_s b_s   \\
    \hline
    0               &     -2k\mu_s a_s   & 2\mu_s\Omega_s \\
    \end{array}
    \right]^{-1}
    \left[
    \begin{array}{cc|c}
    -a_s   &   -k    &   -a_l \\
    -2\mu_s\Omega_s & -2k\mu_s b_s  & -\rho_l\omega^2\\
    \hline
    -2k\mu_s a_s & -2\mu_s\Omega_s & 0 \\
    \end{array}
    \right]

固体-液体界面
----------------

当上层为固体，下层为液体时，推导过程和 `液体-固体界面`_ 完全一致。这里仅给出最终的 R/T 矩阵满足的计算式，读者稍作推导即可验证，

.. math::
    :label:

    \begin{bmatrix}
    \mathbf{T}_D^{1\times2}  & \mathbf{R}_U^{1\times1} \\
    \mathbf{R}_D^{2\times2}  & \mathbf{T}_U^{2\times1} \\
    \end{bmatrix}_{3\times3} = 
    \left[
    \begin{array}{c|cc}
    a_l             &     a_s            &  k             \\
    \rho_l\omega^2  &     2\mu_s\Omega_s &  2k\mu_s b_s   \\
    \hline
    0               &     2k\mu_s a_s    &  2\mu_s\Omega_s \\
    \end{array}
    \right]^{-1}
    \left[
    \begin{array}{cc|c}
    a_s             &   -k            &   a_l \\
    -2\mu_s\Omega_s & 2k\mu_s b_s     &   -\rho_l\omega^2\\
    \hline
    2k\mu_s a_s     & -2\mu_s\Omega_s &   0 \\
    \end{array}
    \right]

之后的操作如增加时间延迟因子，广义 R/T 矩阵递推等不受影响。在程序中为保持 2x2 矩阵，只将以上结果填充在对应位置即可，其余项为0。
