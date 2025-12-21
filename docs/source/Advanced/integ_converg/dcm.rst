:author: 朱邓达
:date: 2025-12-15

直接收敛法（Direct Convergence Method, DCM）
==========================================================

DCM 的想法非常简单（Zhu et al., 2026, BSSA），使用以下一个式子即可说明。
DCM 将原积分拆分为：

.. math::

    \int_0^\infty F(k) J_m(kr) k dk = \int_0^\infty \left[F(k) - F(k_\text{max}) \right] J_m(kr) k dk + F(k_\text{max}) \int_0^\infty J_m(kr) k dk

:math:`k_\text{max}` 即为数值方法的积分上限，而上式第二项积分有如下解析解
，推导见 Zhu et al. (2026) 文章附录，

.. math::

    \int_0^\infty J_m(kr) k dk = 
    \left\{
    \begin{aligned}
    & 0,   & \qquad m=0 \\
    & \dfrac{1}{r^2},  & \qquad m=1 \\
    & \dfrac{2}{r^2},  & \qquad m=2 \\
    \end{aligned}
    \right.

剩下的第一项积分则按照既定的数值方法求解即可。
从以上可见 DCM 几乎不需要额外的复杂计算，且在代码实现上非常简单。

DCM 旨在解决核函数在源点场点深度接近时衰减慢的问题，这与 PTAM 直接作用于振荡的积分结果不同。
以下以半空间模型为例，绘制了场点位于地表时不同震源深度（0 ~ 0.4 km）的 5 Hz 的核函数，
可见浅源的核函数衰减很慢，这是造成积分持续震荡难以收敛的根本原因。

计算 Bash 脚本: :download:`run.sh <run_dcm/run.sh>`
绘图 Python 脚本: :download:`plot_depth_kernel.py <run_dcm/plot_depth_kernel.py>`

.. figure:: run_dcm/deep_shallow_kernel.png
    :align: center

    归一化到 [-1e8, 1e8] 范围内， 使用对称对数轴