:author: 朱邓达
:date: 2026-09-09

平面入射 P-SV 波接收函数的计算公式
========================================

以下介绍 PyGRT 计算的平面入射 P-SV 接收函数所使用的公式。
对于已经熟悉点源计算的读者，可以把两种问题的异同概括为：

* 点源计算对水平波数和柱面谐展开进行积分，而平面波计算固定一个水平波矢，直接指定入射 P 波或 SV 波
* 两者都先由势函数系数得到 :math:`q,w`，再通过已有的分层 R/T 矩阵求解得到自由表面响应

:math:`q/w` 与平面波势函数
----------------------------------------

在 |yao2026p| 中，势函数与位移垂直波函数的关系已经给出。
这里省略点源柱面谐展开的阶数下标 :math:`m`，沿用同一归一化记号：

.. math::
    :label: rcvfn-qw

    q=k\phi+\frac{\partial\psi}{\partial z},
    \qquad
    w=\frac{\partial\phi}{\partial z}+k\psi.

其中 :math:`\phi` 是 P 波势函数，:math:`\Psi` 是双旋度表达式中的原始 SV 标量势。
为了与 P 波势函数采用相同的波数尺度，R/T 矩阵使用缩放后的 SV 势函数
:math:`\psi`，定义为

.. math::

    \psi=k\Psi.

因此，:math:`q,w` 就是底层 R/T 矩阵已经给出的两个响应系数。

平面波固定水平波矢为

.. math::

    \mathbf{k}_h=(k,0),
    \qquad
    \mathcal{E}(x,t)=\exp[\mathrm{i}(kx-\omega t)],
    \qquad
    p=\frac{k}{\omega}.

若在入射层以从竖直方向量起的入射角 :math:`\theta` 指定入射波，则
:math:`p=\sin\theta/V_\alpha`，其中 :math:`V_\alpha` 是对应入射波的速度。
在程序使用射线参数输入时，直接由 :math:`k=\omega p` 得到水平波数。

由势函数得到位移
------------------

位移仍由 P 势函数的一次梯度和 SV 势函数的双旋度组成：

.. math::
    :label: rcvfn-potential

    \mathbf{u}
    =\nabla\phi
    +\nabla\times\left[\nabla\times(\Psi\mathbf{e}_z)\right],
    \qquad
    \psi=k\Psi.

平面波没有 :math:`y` 方向变化，因此

.. math::

    \nabla\times\left[\nabla\times(\Psi\mathbf{e}_z)\right]
    =
    \left(
        \frac{\partial^2\Psi}{\partial x\partial z},
        0,
        -\frac{\partial^2\Psi}{\partial x^2}
    \right).

利用 :math:`\partial_x\exp(\mathrm{i}kx)=\mathrm{i}k\exp(\mathrm{i}kx)`，
水平和垂向位移分量分别为

.. math::
    :label: rcvfn-displacement

    u_R
    &=\frac{\partial\phi}{\partial x}
      +\frac{\partial^2\Psi}{\partial x\partial z}
      =\mathrm{i}\left(k\phi+\frac{\partial\psi}{\partial z}\right)
      =\mathrm{i}q,\\
    u_Z
    &=\frac{\partial\phi}{\partial z}
      -\frac{\partial^2\Psi}{\partial x^2}
      =\frac{\partial\phi}{\partial z}+k\psi
      =w.

所以，平面波计算中真正的位移转换关系是

.. math::

    \boxed{u_R=\mathrm{i}q},
    \qquad
    \boxed{u_Z=w}.

自由表面响应与接收函数
------------------------

完成与点源相同的 R/T 矩阵递推后，在自由表面得到 :math:`q,w`。
记两种入射波对应的自由表面响应为
:math:`(q_P,w_P)` 和 :math:`(q_{SV},w_{SV})`。
直接使用 :eq:`rcvfn-displacement`，即可得到对应分量响应，
再按如下式子即可计算得到接收函数。

P 入射时使用水平分量与垂向分量之比：

.. math::
    :label: rcvfn-ratio-p

    R_P(\omega)
    =\frac{U_R^{(P)}}{U_Z^{(P)}}
    =\frac{\mathrm{i}q_P}{w_P}.

SV 入射时使用垂向分量与水平分量之比：

.. math::
    :label: rcvfn-ratio-sv

    R_{SV}(\omega)
    =\frac{U_Z^{(SV)}}{U_R^{(SV)}}
    =\frac{w_{SV}}{\mathrm{i}q_{SV}}
    =-\frac{\mathrm{i}w_{SV}}{q_{SV}}.

共同的入射振幅在比值中会自动消去，因此计算接收函数时只需使用
R/T 递推得到的 :math:`q,w` 响应。

若还需要输出绝对位移分量，令 :math:`\widetilde V_\alpha` 表示程序在入射层采用的速度因子，
无衰减时它就是 :math:`V_\alpha`。单位物理位移对应的势函数系数为

.. math::

    \Phi_{P,\mathrm{unit}}
    =\frac{\widetilde V_P}{\mathrm{i}\omega},
    \qquad
    \psi_{SV,\mathrm{unit}}
    =\frac{\widetilde V_S}{\omega}
    =\mathrm{i}\frac{\widetilde V_S}{\mathrm{i}\omega}.

因此绝对分量只需在 :eq:`rcvfn-displacement` 的基础上乘以相应的单位入射位移因子；
这不会改变上面的接收函数比值。
