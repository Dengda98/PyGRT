:author: 朱邓达
:date: 2026-08-14

直接收敛法的补充公式
============================================

在应用直接收敛法计算波数积分时，|dcm2025p| 附录给出了基本的校正项公式。
通常将被积函数中的缓慢衰减部分提取出来，
并对剩余部分进行数值积分。设 :math:`F(k)` 是核函数，
:math:`F_\text{max}=F(k_\text{max})`，则最简单的分解为

.. math::

    \int_0^\infty F(k)J_m(kr)k\,\mathrm{d}k
    =
    \int_0^\infty [F(k)-F_\text{max}]J_m(kr)k\,\mathrm{d}k
    +
    F_\text{max}\int_0^\infty J_m(kr)k\,\mathrm{d}k.

本文统一将最后提取部分的解析贡献整体称为DCM 的“校正项”。若该贡献写成
:math:`s \times Q_m`，其中 :math:`s` 表示 :math:`F_\text{max}`、
:math:`F_\text{max}k_\text{max}` 或其他相关系数，:math:`Q_m` 则称为“待求波数积分”，
其具有解析表达式。因此，“校正项”指 :math:`s \times Q_m` 的整体，而不是其中单独的积分 :math:`Q_m`。

不同震源类型的核函数具有不同的高波数衰减形式，因此对不同震源类型设计对应的分解方式能提供更好的效果：

+ **对于单力源，需要提取的部分是** :math:`\dfrac{k_\text{max}}{k} \times F_\text{max}`
+ **对于力偶源，需要提取的部分则是常数** :math:`F_\text{max}`
+ 此外，近场项、位移对 :math:`r` 的偏导以及位移对 :math:`z` 的偏导，**会分别产生不同的波数积分**

本文对这些校正项逐一求解。

以下推导中均假定 :math:`r>0`。下面列出 |dcm2025p| 中使用的相关公式的简化形式：

.. math::

    \int_0^\infty J_\nu(kr)k^\nu\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    \dfrac{1}{r}, &\quad \nu=0,\\
    \dfrac{1}{r^2}, &\quad \nu=1,\\
    \dfrac{3}{r^3}, &\quad \nu=2,
    \end{aligned}
    \right.

.. math::

    \int_0^\infty J_\nu(kr)k^{\nu+1}\,\mathrm{d}k=0,
    \qquad \nu=0,1,2.

因此可得到以下推导中将使用的公式：

.. math::

    \int_0^\infty J_m(kr)\,\mathrm{d}k=\dfrac{1}{r},
    \qquad m=0,1,2,

.. math::

    \int_0^\infty J_m(kr)k\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    0, &\quad m=0,\\
    \dfrac{1}{r^2}, &\quad m=1,\\
    \dfrac{2}{r^2}, &\quad m=2,
    \end{aligned}
    \right.

.. math::

    \int_0^\infty J_m(kr)k^2\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    -\dfrac{1}{r^3}, &\quad m=0,\\
    0, &\quad m=1,\\
    \dfrac{3}{r^3}, &\quad m=2.
    \end{aligned}
    \right.

**以下推导的公式均基于数值积分的方式进行了验证，详见**
:download:`check_DCM_formula.ipynb` （ :doc:`预览 <check_DCM_formula>` ）

波数积分中远场项和近场项的形式
--------------------------------

为说明后续推导中各个波数积分的来源，先写出一个 :math:`m>0` 的典型积分：

.. math::

    \begin{aligned}
    \mathcal I_m
    =
    \int_0^\infty \bigg[
    &q_m(k)J_{m-1}(kr)k
    -m\big(q_m(k)+v_m(k)\big)\frac{J_m(kr)}{r}\\
    &+w_m(k)J_m(kr)k
    -v_m(k)J_{m-1}(kr)k
    \bigg]\,\mathrm{d}k.
    \end{aligned}

其中第一、三、四项属于远场形式，第二项属于近场形式（因为有系数 :math:`1/r`）。对于 :math:`m=0`，
积分中只保留

.. math::

    -\int_0^\infty q_0(k)J_1(kr)k\,\mathrm{d}k
    +
    \int_0^\infty w_0(k)J_0(kr)k\,\mathrm{d}k.

因此，求解 DCM 校正项时，只需分别求解远场的
:math:`F(k)J_m(kr)k` 和近场的 :math:`F(k)J_m(kr)/r`，
以及空间偏导所对应的两个形式。

计算位移
---------

位移积分同时包含远场项和近场项。下面分别按照单力源和力偶源讨论两种情形。

单力源
~~~~~~

远场项
^^^^^^

单力源的高波数核函数按照 :math:`F_\text{max}k_\text{max}/k` 提取，于是

.. math::

    \begin{aligned}
    \int_0^\infty F(k)J_m(kr)k\,\mathrm{d}k
    &=
    \int_0^\infty
    \left[F(k)-F_\text{max}\dfrac{k_\text{max}}{k}\right]
    J_m(kr)k\,\mathrm{d}k\\
    &\quad+
    F_\text{max}k_\text{max}
    \underbrace{\int_0^\infty J_m(kr)\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

其中的待求波数积分在前文已经列出，其结果为

.. math::

    \int_0^\infty J_m(kr)\,\mathrm{d}k=\dfrac{1}{r},
    \qquad m=0,1,2.

近场项
^^^^^^

单力源的高波数核函数按照 :math:`F_\text{max}k_\text{max}/k` 提取，于是

.. math::

    \begin{aligned}
    \dfrac{1}{r}\int_0^\infty F(k)J_m(kr)\,\mathrm{d}k
    &=
    \dfrac{1}{r}\int_0^\infty
    \left[F(k)-F_\text{max}\dfrac{k_\text{max}}{k}\right]
    J_m(kr)\,\mathrm{d}k\\
    &\quad+
    F_\text{max}k_\text{max}
    \underbrace{\dfrac{1}{r}\int_0^\infty
    \dfrac{J_m(kr)}{k}\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

令 :math:`x=kr`，并利用

.. math::

    J_{m-1}(x)+J_{m+1}(x)=\dfrac{2m}{x}J_m(x),
    \qquad
    \int_0^\infty J_m(x)\,\mathrm{d}x=1,
    \qquad m=1,2,

可得

.. math::

    \int_0^\infty \dfrac{J_m(x)}{x}\,\mathrm{d}x
    =\dfrac{1}{2m}\int_0^\infty
    [J_{m-1}(x)+J_{m+1}(x)]\,\mathrm{d}x
    =\dfrac{1}{m},
    \qquad m=1,2.

因此单力源近场项的待求波数积分为

.. math::

    \dfrac{1}{r}\int_0^\infty \dfrac{J_m(kr)}{k}\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    \dfrac{1}{r}, &\quad m=1,\\
    \dfrac{1}{2r}, &\quad m=2.
    \end{aligned}
    \right.

力偶源
~~~~~~

远场项
^^^^^^

力偶源的高波数核函数提取为常数，于是

.. math::

    \begin{aligned}
    \int_0^\infty F(k)J_m(kr)k\,\mathrm{d}k
    &=
    \int_0^\infty [F(k)-F_\text{max}]J_m(kr)k\,\mathrm{d}k\\
    &\quad+
    F_\text{max}
    \underbrace{\int_0^\infty J_m(kr)k\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

其中的待求波数积分在前文已经给出，可直接使用相应结果：

.. math::

    \int_0^\infty J_m(kr)k\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    0, &\quad m=0,\\
    \dfrac{1}{r^2}, &\quad m=1,\\
    \dfrac{2}{r^2}, &\quad m=2.
    \end{aligned}
    \right.

近场项
^^^^^^

力偶源的高波数核函数提取为常数，于是

.. math::

    \begin{aligned}
    \dfrac{1}{r}\int_0^\infty F(k)J_m(kr)\,\mathrm{d}k
    &=
    \dfrac{1}{r}\int_0^\infty
    [F(k)-F_\text{max}]J_m(kr)\,\mathrm{d}k\\
    &\quad+
    F_\text{max}
    \underbrace{\dfrac{1}{r}\int_0^\infty J_m(kr)\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

其中的积分部分与前面列出的波数积分结果一致，近场项额外带有 :math:`1/r` 因子，因此

.. math::

    \dfrac{1}{r}\int_0^\infty J_m(kr)\,\mathrm{d}k
    =\dfrac{1}{r^2},
    \qquad m=1,2.

计算位移对 :math:`r` 的偏导
------------------------------------

位移对 :math:`r` 求偏导时，远场项和近场项分别变为

.. math::

    \int_0^\infty F(k)J_m'(kr)k^2\,\mathrm{d}k,
    \qquad
    \int_0^\infty F(k)D_m(k,r)\,\mathrm{d}k,

其中

.. math::

    D_m(k,r)
    =
    \dfrac{\mathrm{d}}{\mathrm{d}r}
    \left[\dfrac{1}{r}J_m(kr)\right]
    =
    \dfrac{k}{r}J_m'(kr)-\dfrac{1}{r^2}J_m(kr).

单力源
~~~~~~

远场项
^^^^^^

按照 :math:`F_\text{max}k_\text{max}/k` 提取，有

.. math::

    \begin{aligned}
    \int_0^\infty F(k)J_m'(kr)k^2\,\mathrm{d}k
    &=
    \int_0^\infty
    \left[F(k)-F_\text{max}\dfrac{k_\text{max}}{k}\right]
    J_m'(kr)k^2\,\mathrm{d}k\\
    &\quad+
    F_\text{max}k_\text{max}
    \underbrace{\int_0^\infty J_m'(kr)k\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

利用递推关系

.. math::

    J_m'(x)=J_{m-1}(x)-\dfrac{m}{x}J_m(x),

分别得到

.. math::

    \begin{aligned}
    \int_0^\infty J_0'(kr)k\,\mathrm{d}k
    &=-\int_0^\infty J_1(kr)k\,\mathrm{d}k=-\dfrac{1}{r^2},\\
    \int_0^\infty J_1'(kr)k\,\mathrm{d}k
    &=\dfrac{1}{r}\int_0^\infty J_1(kr)\,\mathrm{d}k
      -\int_0^\infty J_2(kr)k\,\mathrm{d}k\\
    &=\dfrac{1}{r^2}-\dfrac{2}{r^2}=-\dfrac{1}{r^2},\\
    \int_0^\infty J_2'(kr)k\,\mathrm{d}k
    &=\int_0^\infty J_1(kr)k\,\mathrm{d}k
      -\dfrac{2}{r}\int_0^\infty J_2(kr)\,\mathrm{d}k\\
    &=\dfrac{1}{r^2}-\dfrac{2}{r^2}=-\dfrac{1}{r^2}.
    \end{aligned}

所以单力源远场项的待求波数积分为

.. math::

    \int_0^\infty J_m'(kr)k\,\mathrm{d}k
    =-\dfrac{1}{r^2},
    \qquad m=0,1,2.

近场项
^^^^^^

近场项的提取形式与位移的近场项相同，因此

.. math::

    \begin{aligned}
    \int_0^\infty F(k)D_m(k,r)\,\mathrm{d}k
    &=
    \int_0^\infty
    \left[F(k)-F_\text{max}\dfrac{k_\text{max}}{k}\right]
    D_m(k,r)\,\mathrm{d}k\\
    &\quad+
    F_\text{max}k_\text{max}
    \underbrace{\int_0^\infty \dfrac{D_m(k,r)}{k}\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

由 :math:`D_m(k,r)` 的定义，

.. math::

    \int_0^\infty \dfrac{D_m(k,r)}{k}\,\mathrm{d}k
    =
    \dfrac{1}{r}\int_0^\infty J_m'(kr)\,\mathrm{d}k
    -\dfrac{1}{r^2}\int_0^\infty \dfrac{J_m(kr)}{k}\,\mathrm{d}k.

其中第一项为零，第二项使用前面位移的单力源近场结果，得到

.. math::

    \int_0^\infty \dfrac{D_m(k,r)}{k}\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    -\dfrac{1}{r^2}, &\quad m=1,\\
    -\dfrac{1}{2r^2}, &\quad m=2.
    \end{aligned}
    \right.

力偶源
~~~~~~

远场项
^^^^^^

力偶源的分解为

.. math::

    \begin{aligned}
    \int_0^\infty F(k)J_m'(kr)k^2\,\mathrm{d}k
    &=
    \int_0^\infty [F(k)-F_\text{max}]
    J_m'(kr)k^2\,\mathrm{d}k\\
    &\quad+
    F_\text{max}
    \underbrace{\int_0^\infty J_m'(kr)k^2\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

由同一递推关系

.. math::

    \begin{aligned}
    \int_0^\infty J_0'(kr)k^2\,\mathrm{d}k
    &=-\int_0^\infty J_1(kr)k^2\,\mathrm{d}k=0,\\
    \int_0^\infty J_1'(kr)k^2\,\mathrm{d}k
    &=\dfrac{1}{r}\int_0^\infty J_1(kr)k\,\mathrm{d}k
      -\int_0^\infty J_2(kr)k^2\,\mathrm{d}k\\
    &=\dfrac{1}{r^3}-\dfrac{3}{r^3}=-\dfrac{2}{r^3},\\
    \int_0^\infty J_2'(kr)k^2\,\mathrm{d}k
    &=\int_0^\infty J_1(kr)k^2\,\mathrm{d}k
      -\dfrac{2}{r}\int_0^\infty J_2(kr)k\,\mathrm{d}k\\
    &=0-\dfrac{4}{r^3}=-\dfrac{4}{r^3}.
    \end{aligned}

所以力偶源远场项的待求波数积分为

.. math::

    \int_0^\infty J_m'(kr)k^2\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    0, &\quad m=0,\\
    -\dfrac{2}{r^3}, &\quad m=1,\\
    -\dfrac{4}{r^3}, &\quad m=2.
    \end{aligned}
    \right.

近场项
^^^^^^

近场项的提取形式与位移的近场项相同，但这里的微分对象为 :math:`D_m(k,r)`，因此

.. math::

    \begin{aligned}
    \int_0^\infty F(k)D_m(k,r)\,\mathrm{d}k
    &=
    \int_0^\infty [F(k)-F_\text{max}]D_m(k,r)\,\mathrm{d}k\\
    &\quad+
    F_\text{max}\underbrace{\int_0^\infty D_m(k,r)\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

由 :math:`D_m(k,r)` 的定义，并沿用前面已经得到的波数积分结果，可得

.. math::

    \int_0^\infty D_m(k,r)\,\mathrm{d}k
    =
    \dfrac{1}{r}\int_0^\infty J_m'(kr)k\,\mathrm{d}k
    -\dfrac{1}{r^2}\int_0^\infty J_m(kr)\,\mathrm{d}k
    =-\dfrac{2}{r^3},
    \qquad m=1,2.

计算位移对 :math:`z` 的偏导
------------------------------------

位移对 :math:`z` 的偏导作用在核函数上，因此 Bessel 函数部分不发生微分，
远场项和近场项仍分别保持前面给出的形式。

单力源
~~~~~~

远场项
^^^^^^

单力源的 :math:`z` 偏导核函数趋于常数，因此

.. math::

    \begin{aligned}
    \int_0^\infty F_z(k)J_m(kr)k\,\mathrm{d}k
    &=
    \int_0^\infty [F_z(k)-F_{z,\text{max}}]
    J_m(kr)k\,\mathrm{d}k\\
    &\quad+
    F_{z,\text{max}}
    \underbrace{\int_0^\infty J_m(kr)k\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

这里的待求波数积分与位移部分中力偶源的远场项完全一致，可直接沿用其结果：

.. math::

    \int_0^\infty J_m(kr)k\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    0, &\quad m=0,\\
    \dfrac{1}{r^2}, &\quad m=1,\\
    \dfrac{2}{r^2}, &\quad m=2.
    \end{aligned}
    \right.

近场项
^^^^^^

单力源近场项的分解为

.. math::

    \begin{aligned}
    \dfrac{1}{r}\int_0^\infty F_z(k)J_m(kr)\,\mathrm{d}k
    &=
    \dfrac{1}{r}\int_0^\infty
    [F_z(k)-F_{z,\text{max}}]J_m(kr)\,\mathrm{d}k\\
    &\quad+
    F_{z,\text{max}}
    \underbrace{\dfrac{1}{r}\int_0^\infty J_m(kr)\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

这里的待求波数积分与位移部分中力偶源的近场项一致，可直接沿用其结果：

.. math::

    \dfrac{1}{r}\int_0^\infty J_m(kr)\,\mathrm{d}k
    =\dfrac{1}{r^2},
    \qquad m=1,2.

力偶源
~~~~~~

远场项
^^^^^^

力偶源的 :math:`z` 偏导核函数按照
:math:`F_{z,\text{max}}k/k_\text{max}` 提取，因此

.. math::

    \begin{aligned}
    \int_0^\infty F_z(k)J_m(kr)k\,\mathrm{d}k
    &=
    \int_0^\infty
    k\left[\dfrac{F_z(k)}{k}
    -\dfrac{F_{z,\text{max}}}{k_\text{max}}\right]
    J_m(kr)k\,\mathrm{d}k\\
    &\quad+
    \dfrac{F_{z,\text{max}}}{k_\text{max}}
    \underbrace{\int_0^\infty J_m(kr)k^2\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

这里的待求波数积分就是前面列出的波数积分结果：

.. math::

    \int_0^\infty J_m(kr)k^2\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    -\dfrac{1}{r^3}, &\quad m=0,\\
    0, &\quad m=1,\\
    \dfrac{3}{r^3}, &\quad m=2.
    \end{aligned}
    \right.

近场项
^^^^^^

力偶源近场项的分解为

.. math::

    \begin{aligned}
    \dfrac{1}{r}\int_0^\infty F_z(k)J_m(kr)\,\mathrm{d}k
    &=
    \dfrac{1}{r}\int_0^\infty
    k\left[\dfrac{F_z(k)}{k}
    -\dfrac{F_{z,\text{max}}}{k_\text{max}}\right]
    J_m(kr)\,\mathrm{d}k\\
    &\quad+
    \dfrac{F_{z,\text{max}}}{k_\text{max}}
    \underbrace{\dfrac{1}{r}\int_0^\infty J_m(kr)k\,\mathrm{d}k}_{\text{待求波数积分}}.
    \end{aligned}

其中的积分部分与位移部分中力偶源远场项的待求波数积分相同，近场项额外带有 :math:`1/r` 因子，因此

.. math::

    \dfrac{1}{r}\int_0^\infty J_m(kr)k\,\mathrm{d}k
    =
    \left\{
    \begin{aligned}
    \dfrac{1}{r^3}, &\quad m=1,\\
    \dfrac{2}{r^3}, &\quad m=2.
    \end{aligned}
    \right.
