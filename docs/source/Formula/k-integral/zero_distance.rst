:author: 朱邓达
:date: 2026-08-15

零震中距的处理
============================================

当场点与震源的水平距离为零时，即震中距 :math:`r=0`，此时，Bessel 函数值为常数，
且含有 :math:`1/r` 的波数积分表达式不能直接计算，需要先取 Bessel 函数的极限。
在程序计算中，:math:`r=0` 作为独立情形处理。

Bessel 函数在 :math:`r=0` 附近的展开式为

.. math::

    J_0(kr)=1-\dfrac{k^2r^2}{4}+O(r^4),
    \qquad
    J_1(kr)=\dfrac{kr}{2}+O(r^3),
    \qquad
    J_2(kr)=\dfrac{k^2r^2}{8}+O(r^4).

因此

.. math::

    J_0(0)=1,
    \qquad
    J_1(0)=J_2(0)=0.

对于近场表达式中出现的 :math:`J_m(kr)/r`，有

.. math::

    \lim_{r\to0}\dfrac{J_1(kr)}{r}=\dfrac{k}{2},
    \qquad
    \lim_{r\to0}\dfrac{J_2(kr)}{r}=0.

对上述近场表达式继续对 :math:`r` 求导，可得

.. math::

    \lim_{r\to0}
    \dfrac{\mathrm{d}}{\mathrm{d}r}
    \left[\dfrac{J_1(kr)}{r}\right]=0,
    \qquad
    \lim_{r\to0}
    \dfrac{\mathrm{d}}{\mathrm{d}r}
    \left[\dfrac{J_2(kr)}{r}\right]=\dfrac{k^2}{8}.

波数积分在 :math:`r=0` 时的形式
------------------------------------

设 :math:`F(k)` 为任意核函数。以下各式中的 :math:`I(0)` 均表示
:math:`r\to0` 时的极限值。对于远场形式的波数积分

.. math::

    I_m^{\mathrm{far}}(r)
    =
    \int_0^\infty F(k)J_m(kr)k\,\mathrm{d}k,

利用 :math:`J_0(0)=1` 以及 :math:`J_1(0)=J_2(0)=0`，可得

.. math::

    I_m^{\mathrm{far}}(0)
    =
    \left\{
    \begin{aligned}
    \int_0^\infty F(k)k\,\mathrm{d}k, &\quad m=0,\\
    0, &\quad m=1,2.
    \end{aligned}
    \right.

对于近场形式的波数积分，:math:`m=1,2` 时有

.. math::

    I_m^{\mathrm{near}}(r)
    =
    \dfrac{1}{r}\int_0^\infty F(k)J_m(kr)\,\mathrm{d}k,

因此

.. math::

    I_m^{\mathrm{near}}(0)
    =
    \left\{
    \begin{aligned}
    \dfrac{1}{2}\int_0^\infty F(k)k\,\mathrm{d}k, &\quad m=1,\\
    0, &\quad m=2.
    \end{aligned}
    \right.

对远场项求 :math:`r` 偏导时，积分形式为

.. math::

    \dfrac{\partial I_m^{\mathrm{far}}}{\partial r}
    =
    \int_0^\infty F(k)J_m'(kr)k^2\,\mathrm{d}k.

由 :math:`J_0'(0)=J_2'(0)=0` 和 :math:`J_1'(0)=1/2`，得到

.. math::

    \left.\dfrac{\partial I_m^{\mathrm{far}}}{\partial r}\right|_{r=0}
    =
    \left\{
    \begin{aligned}
    0, &\quad m=0,2,\\
    \dfrac{1}{2}\int_0^\infty F(k)k^2\,\mathrm{d}k, &\quad m=1.
    \end{aligned}
    \right.

对近场项求 :math:`r` 偏导时，利用前面的极限可得

.. math::

    \left.\dfrac{\partial I_m^{\mathrm{near}}}{\partial r}\right|_{r=0}
    =
    \left\{
    \begin{aligned}
    0, &\quad m=1,\\
    \dfrac{1}{8}\int_0^\infty F(k)k^2\,\mathrm{d}k, &\quad m=2.
    \end{aligned}
    \right.

如果空间偏导作用在核函数上，例如将 :math:`F(k)` 替换为 :math:`F_z(k)`，
上述远场和近场结果的 Bessel 函数部分不变，只需将对应的核函数代入这些普通波数积分即可。
这些表达式成立的前提是取极限后得到的 :math:`k` 积分本身收敛。
