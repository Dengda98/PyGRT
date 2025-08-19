反射透射系数矩阵的另一种推导
==============================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

在 :ref:`《理论地震图及其应用（初稿）》 <yao_init_manuscripts>` 中给出了反射透射系数矩阵（下称 R/T 矩阵）的推导，其中推导的基础是先计算出固体-固体界面两侧的垂直波函数转换矩阵 :math:`\mathbf{Q}(j^-,j^+) = \mathbf{D}_{j-1}^{-1} \mathbf{D}_{j}` ，但这种推导不太适用于后续对于液体-固体界面的 R/T 矩阵的推导，故这里给出另一种推导方式，本质还是矩阵的各种变换。

原推导过程涉及多种矢量符号，比较容易产生混淆，以下尽量将矩阵、矢量展开。

以下推导以动态 P-SV 情况为例。

最终每一项系数的详细表达式我使用 Python 库 `SymPy <https://www.sympy.org/>`_ 辅助推导，包括与原推导结果的对比验证，可在这里下载 :download:`RT_formula.ipynb` （ :doc:`预览 <RT_formula>` ）。

----------------------------------------


界面两侧的垂直波函数关系
-----------------------------

从波动方程出发，使用柱面谐矢量(vector cylindrical harmonics)表示位移和应力，则某层内的位移垂直波函数与位移应力矢量之间的关系可以表示为（略去角标 :math:`j`）

.. math::
    :label:

    \begin{bmatrix}
    q_m \\
    w_m \\
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

其中左侧为位移应力矢量，右侧的 4x4 矩阵 我们称为 :math:`D` 矩阵。为简单起见，我们假设界面位于第1/2层之间，根据应力位移连续的边界条件，有

.. math::
    :label: layer12

    & \begin{bmatrix}
    k     &   b_1   &   k    &  -b_1 \\
    a_1   &   k     &  -a_1  &  k  \\
    2\mu_1\Omega_1 & 2k\mu_1 b_1 & 2\mu_1\Omega_1 & -2k\mu_1 b_1 \\
    2k\mu_1 a_1 & 2\mu_1\Omega_1 & -2k\mu_1 a_1 & 2\mu_1\Omega_1 \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^- (z_1^-) \\
    \psi_m^- (z_1^-) \\
    \phi_m^+ (z_1^-) \\
    \psi_m^+ (z_1^-) \\
    \end{bmatrix} \\
    = & \begin{bmatrix}
    k     &   b_2   &   k  &  -b_2 \\
    a_2   &   k     &  -a_2  &  k  \\
    2\mu_2\Omega_2 & 2k\mu_2 b_2 & 2\mu_2\Omega_2 & -2k\mu_2 b_2 \\
    2k\mu_2 a_2 & 2\mu_2\Omega_2 & -2k\mu_2 a_2 & 2\mu_2\Omega_2 \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^- (z_1^+) \\
    \psi_m^- (z_1^+) \\
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{bmatrix}

界面上的 R/T 矩阵
--------------------------

原推导中对 :eq:`layer12` 式做了矩阵求逆以获得两侧垂直波函数的转换矩阵。这里我们先保留该形式，仍然根据入射方向分两种情况讨论，最终同样可求出 R/T 矩阵。

波从上向下入射
~~~~~~~~~~~~~~~~

此时下层没有向上传播的入射波，即 :math:`[\phi_m^- (z_1^+), \psi_m^- (z_1^+)]^T = \mathbf{0}` ，:eq:`layer12` 式变为

.. math::
    :label:

    \begin{bmatrix}
    k     &   b_1   &   k    &  -b_1 \\
    a_1   &   k     &  -a_1  &  k  \\
    2\mu_1\Omega_1 & 2k\mu_1 b_1 & 2\mu_1\Omega_1 & -2k\mu_1 b_1 \\
    2k\mu_1 a_1 & 2\mu_1\Omega_1 & -2k\mu_1 a_1 & 2\mu_1\Omega_1 \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^- (z_1^-) \\
    \psi_m^- (z_1^-) \\
    \bbox[yellow] {\phi_m^+ (z_1^-)} \\
    \bbox[yellow] {\psi_m^+ (z_1^-)} \\
    \end{bmatrix} = \begin{bmatrix}
    k  &  -b_2 \\
    -a_2  &  k  \\
    2\mu_2\Omega_2 & -2k\mu_2 b_2 \\
    -2k\mu_2 a_2 & 2\mu_2\Omega_2 \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{bmatrix}

其中高亮部分的垂直波函数为当前情况的“已知项”，通过移项+矩阵重排的方式可得到

.. math::
    :label: U2D

    \begin{bmatrix}
    -k     &   -b_1   &   k    &  -b_2 \\
    -a_1   &   -k     &  -a_2  &  k  \\
    -2\mu_1\Omega_1 & -2k\mu_1 b_1 & 2\mu_2\Omega_2 & -2k\mu_2 b_2 \\
    -2k\mu_1 a_1 & -2\mu_1\Omega_1 & -2k\mu_2 a_2 & 2\mu_2\Omega_2 \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^- (z_1^-) \\
    \psi_m^- (z_1^-) \\
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{bmatrix} = 
    \begin{bmatrix}
    k  &  -b_1 \\
    -a_1  &  k  \\
    2\mu_1\Omega_1 & -2k\mu_1 b_1 \\
    -2k\mu_1 a_1 & 2\mu_1\Omega_1 \\
    \end{bmatrix} 
    \begin{bmatrix}
    \bbox[yellow] {\phi_m^+ (z_1^-)} \\
    \bbox[yellow] {\psi_m^+ (z_1^-)} \\
    \end{bmatrix} 

其中等号左边矩阵前两列的负号由移项产生，此时左边的垂直波函数矢量（作为未知量）已经变成两层的混合版本，适定方程可简单使用逆矩阵求解，得到

.. math::
    :label:

    \begin{bmatrix}
    \phi_m^- (z_1^-) \\
    \psi_m^- (z_1^-) \\
    \end{bmatrix} = 
    \mathbf{R}_D^{2\times2}
    \begin{bmatrix}
    \bbox[yellow] {\phi_m^+ (z_1^-)} \\
    \bbox[yellow] {\psi_m^+ (z_1^-)} \\
    \end{bmatrix} 

    \begin{bmatrix}
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{bmatrix} = 
    \mathbf{T}_D^{2\times2}
    \begin{bmatrix}
    \bbox[yellow] {\phi_m^+ (z_1^-)} \\
    \bbox[yellow] {\psi_m^+ (z_1^-)} \\
    \end{bmatrix} 

波从下向上入射
~~~~~~~~~~~~~~~~

此时上层没有向下传播的入射波，即 :math:`[\phi_m^+ (z_1^-), \psi_m^+ (z_1^-)]^T = \mathbf{0}` ，:eq:`layer12` 式变为

.. math::
    :label:

    \begin{bmatrix}
    k     &   b_1   \\
    a_1   &   k     \\
    2\mu_1\Omega_1 & 2k\mu_1 b_1  \\
    2k\mu_1 a_1 & 2\mu_1\Omega_1  \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^- (z_1^-) \\
    \psi_m^- (z_1^-) \\
    \end{bmatrix} = \begin{bmatrix}
    k     &   b_2   &   k  &  -b_2 \\
    a_2   &   k     &  -a_2  &  k  \\
    2\mu_2\Omega_2 & 2k\mu_2 b_2 & 2\mu_2\Omega_2 & -2k\mu_2 b_2 \\
    2k\mu_2 a_2 & 2\mu_2\Omega_2 & -2k\mu_2 a_2 & 2\mu_2\Omega_2 \\
    \end{bmatrix} 
    \begin{bmatrix}
    \bbox[yellow] {\phi_m^- (z_1^+)} \\
    \bbox[yellow] {\psi_m^- (z_1^+)} \\
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{bmatrix}

其中高亮部分的垂直波函数同样为当前情况的“已知项”，为保持与 :eq:`U2D` 式的形式匹配，通过类似的移项+矩阵重排的方式可得到

.. math::
    :label: D2U

    \begin{bmatrix}
    -k     &   -b_1   &   k    &  -b_2 \\
    -a_1   &   -k     &  -a_2  &  k  \\
    -2\mu_1\Omega_1 & -2k\mu_1 b_1 & 2\mu_2\Omega_2 & -2k\mu_2 b_2 \\
    -2k\mu_1 a_1 & -2\mu_1\Omega_1 & -2k\mu_2 a_2 & 2\mu_2\Omega_2 \\
    \end{bmatrix} 
    \begin{bmatrix}
    \phi_m^- (z_1^-) \\
    \psi_m^- (z_1^-) \\
    \phi_m^+ (z_1^+) \\
    \psi_m^+ (z_1^+) \\
    \end{bmatrix} = 
    \begin{bmatrix}
    -k  &  -b_2 \\
    -a_2  &  -k  \\
    -2\mu_2\Omega_2 & -2k\mu_2 b_2 \\
    -2k\mu_2 a_2 & -2\mu_2\Omega_2 \\
    \end{bmatrix} 
    \begin{bmatrix}
    \bbox[yellow] {\phi_m^- (z_1^+)} \\
    \bbox[yellow] {\psi_m^- (z_1^+)} \\
    \end{bmatrix} 

矩阵中的负号由移项产生，等号左边形式与 :eq:`U2D` 式完全一致。同样该适定方程可简单使用逆矩阵求解，得到

.. math::
    :label:

    \begin{bmatrix}
    \phi_m^- (z_1^-) \\
    \psi_m^- (z_1^-) \\
    \end{bmatrix} = 
    \mathbf{T}_U^{2\times2}
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

:eq:`U2D` 式和 :eq:`D2U` 式可合并，一并使用逆矩阵求得最终界面上的 R/T 矩阵，

.. math::
    :label:

    & \begin{bmatrix}
    \mathbf{T}_U^{2\times2}  & \mathbf{R}_D^{2\times2} \\
    \mathbf{R}_U^{2\times2}  & \mathbf{T}_D^{2\times2} \\
    \end{bmatrix} \\
    = & 
    \begin{bmatrix}
    -k     &   -b_1   &   k    &  -b_2 \\
    -a_1   &   -k     &  -a_2  &  k  \\
    -2\mu_1\Omega_1 & -2k\mu_1 b_1 & 2\mu_2\Omega_2 & -2k\mu_2 b_2 \\
    -2k\mu_1 a_1 & -2\mu_1\Omega_1 & -2k\mu_2 a_2 & 2\mu_2\Omega_2 \\
    \end{bmatrix}^{-1}
    \begin{bmatrix}
    -k     &   -b_2   &   k    &  -b_1 \\
    -a_2   &   -k     &  -a_1  &  k  \\
    -2\mu_2\Omega_2 & -2k\mu_2 b_2 & 2\mu_1\Omega_1 & -2k\mu_1 b_1 \\
    -2k\mu_2 a_2 & -2\mu_2\Omega_2 & -2k\mu_1 a_1 & 2\mu_1\Omega_1 \\
    \end{bmatrix}

之后的运算如增加时间延迟因子，广义 R/T 矩阵递推等不受影响。