.. _static_uniform:

无限均匀介质中集中脉冲力产生的静态位移场
===========================================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

这是一个经典问题，通常使用势函数分解（例如动态位移场的求解）。这里提供另一种解法，空间域到波矢域的积分变换法。结果表达式将用于确定震源系数。

静态弹性波方程
------------------------------
在直角坐标系下，研究位于原点处沿单位矢量 :math:`\mathbf{n}` 作用的一个集中脉冲力 :math:`\delta(\mathbf{r})`，在 :math:`\mathbf{r}=(x_1, x_2, x_3)` 处产生的静态位移场 :math:`\mathbf{u}=(u_1, u_2, u_3)` 。此时弹性波方程表示为

.. math:: 
    :label:

    (\lambda+2\mu) \nabla(\nabla \cdot \mathbf{u}) - \mu \nabla \times (\nabla \times \mathbf{u}) = - \delta(\mathbf{r}) \mathbf{n}

利用 :math:`\nabla \times (\nabla \times A) = \nabla(\nabla \cdot A) - \nabla^2 A` 展开双旋度，

.. math:: 
    :label:

    (\lambda+\mu) \nabla(\nabla \cdot \mathbf{u}) + \mu \nabla^2 \mathbf{u} = - \delta(\mathbf{r}) \mathbf{n}

不妨假设 :math:`\mathbf{n}` 沿着某个基矢量 :math:`\mathbf{e_j}`，并引入格林函数 :math:`\mathbf{G}_j=(G_{1j}, G_{2j}, G_{3j})`，

.. math:: 
    :label:

    (\lambda+\mu) \nabla(\nabla \cdot \mathbf{G}_j) + \mu \nabla^2 \mathbf{G}_j = - \delta(\mathbf{r}) \mathbf{e}_j


从三维空间域变换到三维波矢域，利用微分性质 :math:`\nabla f(\mathbf{r}) \leftrightarrow i\mathbf{k} F(\mathbf{k})`，

.. math:: 
    :label:

    (\lambda+\mu) \mathbf{k}(\mathbf{k} \cdot \mathbf{\tilde{G}}_j) + \mu k^2 \mathbf{\tilde{G}}_j = \mathbf{e}_j

其中 :math:`\mathbf{\tilde{G}}_j` 为 :math:`\mathbf{G}_j` 的Fourier变换，:math:`\mathbf{k}=(k_1, k_2, k_3)` 为波矢量，:math:`k=|\mathbf{k}|`。移项整理得到 :math:`\mathbf{\tilde{G}}_j` 的表达式，

.. math:: 
    :label:

    \mathbf{\tilde{G}}_j = \dfrac{1}{\mu k^2} \left[ \mathbf{I} + (\lambda+\mu) \dfrac{\mathbf{kk^T}}{k^2} \right]^{-1} \mathbf{e}_j

使用Sherman-Morrison公式 [#]_ 计算上式逆矩阵，

.. math:: 
    :label:

    (\mathbf{A} + \mathbf{bc^T}) = \mathbf{A}^{-1} - \dfrac{\mathbf{A}^{-1}\mathbf{bc^T}\mathbf{A}^{-1}}{1 + \mathbf{c^T}\mathbf{A}^{-1}\mathbf{b}}

得到 :math:`\mathbf{\tilde{G}}_j` 的表达式，

.. math:: 
    :label: static_Gk

    \mathbf{\tilde{G}}_j = \dfrac{1}{\mu k^2} \left[ \mathbf{I} - \dfrac{\lambda+\mu}{\lambda+2\mu} \dfrac{\mathbf{kk^T}}{k^2} \right] \mathbf{e}_j
    = \dfrac{1}{\mu k^2} \left[ \mathbf{I} - \dfrac{1}{2(1-\nu)} \dfrac{\mathbf{kk^T}}{k^2} \right] \mathbf{e}_j

其中 :math:`\nu` 为泊松比，:math:`\nu=\dfrac{\lambda}{2(\lambda+\mu)}`。



Fourier变换对
------------------------------
对于三维空间中的函数 :math:`f(\mathbf{r})`，其Fourier变换 [#]_ 定义为

.. math:: 
    :label: 

    F(\mathbf{k}) = \mathcal{F} \left[ f(\mathbf{r}) \right] = \int_{-\infty}^{+\infty} \int_{-\infty}^{+\infty} \int_{-\infty}^{+\infty} f(\mathbf{r}) e^{-i \mathbf{k}\cdot \mathbf{r}} d^3 \mathbf{r}


对应的逆变换为

.. math:: 
    :label: 

    f(\mathbf{r}) = \mathcal{F}^{-1} \left[ F(\mathbf{k}) \right] = \dfrac{1}{(2\pi)^3} \int_{-\infty}^{+\infty} \int_{-\infty}^{+\infty} \int_{-\infty}^{+\infty} F(\mathbf{k}) e^{i \mathbf{k}\cdot \mathbf{r}} d^3 \mathbf{k}


为了将 :eq:`static_Gk` 式变换回空间域，可逐项做逆变换，再合并得到最终解。 不过注意到 :eq:`static_Gk` 式中的特殊形式， **可基于微分关系快速获得对应项的变换对** 。

已知 :math:`1/r` 满足以下关系，

.. math:: 
    :label:

    \nabla^2 \left( \dfrac{1}{r} \right) = - 4\pi \delta(\mathbf{r})

其中 :math:`r=|\mathbf{r}|` 。变换到波矢域，得到变换对，

.. math:: 
    :label: 1k2_inv_final

    \mathcal{F} \left[ \dfrac{1}{4\pi r} \right] &= \dfrac{1}{k^2} \\
    \mathcal{F}^{-1} \left[ \dfrac{1}{k^2}  \right] &= \dfrac{1}{4\pi r}

而 :math:`1/r` 还满足以下变换关系，

.. math:: 
    :label:

    \mathcal{F} \left[ \nabla \nabla \left( \dfrac{1}{4\pi r} \right) \right] = \dfrac{\mathbf{kk^T}}{k^2}


因此 :eq:`static_Gk` 式中涉及对 :math:`\mathbf{kk^T}/k^4` 的逆变换为

.. math:: 
    :label: kkTk4_inv

    \mathcal{F}^{-1} \left[ \dfrac{1}{k^2} \cdot \dfrac{\mathbf{kk^T}}{k^2} \right] = 
    \dfrac{1}{4\pi r} * \nabla \nabla \left( \dfrac{1}{4\pi r} \right) = 
    \dfrac{1}{(4\pi)^2} \cdot \nabla \nabla \left( \dfrac{1}{r} * \dfrac{1}{r} \right) 

同样再次使用积分变换求解卷积式 :math:`\left( \dfrac{1}{r} * \dfrac{1}{r} \right)`，其正变换为

.. math:: 
    :label: 1r1r_conv

    \dfrac{1}{r} * \dfrac{1}{r} = \mathcal{F}^{-1} \left[ \mathcal{F} \left[ \dfrac{1}{r} * \dfrac{1}{r} \right] \right] = 
    \mathcal{F}^{-1} \left[ \mathcal{F} \left[ \dfrac{1}{r} \right] \cdot \mathcal{F} \left[ \dfrac{1}{r} \right] \right] 
    = \mathcal{F}^{-1} \left[ \dfrac{(4\pi)^2}{k^4} \right] 

注意到 :math:`r` 满足以下关系，

.. math:: 
    :label:

    \nabla^2 \left( \nabla^2 r \right) = \nabla^2 \left( \dfrac{2}{r} \right) = 8\pi \delta(\mathbf{r})

变换到波矢域，得到变换对，

.. math:: 
    :label:

    \mathcal{F} \left[ \dfrac{r}{8\pi} \right] &= \dfrac{1}{k^4} \\
    \mathcal{F}^{-1} \left[ \dfrac{1}{k^4} \right] &= \dfrac{r}{8\pi}

上式代入到 :eq:`1r1r_conv` 式，得到

.. math:: 
    :label:

    \dfrac{1}{r} * \dfrac{1}{r} = 2\pi r 

上式代入到 :eq:`kkTk4_inv` 式，得到

.. math:: 
    :label: kkTk4_inv_final

    \mathcal{F}^{-1} \left[ \dfrac{1}{k^2} \cdot \dfrac{\mathbf{kk^T}}{k^2} \right] = 
    \dfrac{1}{(4\pi)^2} \cdot \nabla \nabla \left( 2\pi r \right) =
    \dfrac{1}{8\pi} \left( \dfrac{\mathbf{I}}{r} - \dfrac{\mathbf{rr^T}}{r^3} \right)


静态位移场的最终表达式
-----------------------------------


将 :eq:`1k2_inv_final` 式，:eq:`kkTk4_inv_final` 式代入 :eq:`static_Gk` 式，作逆变换，整理，

.. math:: 
    :label: 

    \mathbf{G}_j &= \dfrac{1}{\mu} \left[ \dfrac{1}{4\pi r} - \dfrac{1}{2(1-\nu)} \dfrac{1}{8\pi} \left( \dfrac{\mathbf{I}}{r} - \dfrac{\mathbf{rr^T}}{r^3} \right) \right]  \mathbf{e}_j \\
    &= \dfrac{1}{16 \pi \mu r} \left[ 4 \mathbf{I} - \dfrac{1}{1-\nu} \mathbf{I} + \dfrac{1}{1-\nu} \dfrac{\mathbf{rr^T}}{r^2}  \right]  \mathbf{e}_j \\
    &= \dfrac{1}{16 \pi \mu (1-\nu) r}  \left[ (3-4\nu)\mathbf{I} + \dfrac{\mathbf{rr^T}}{r^2} \right]  \mathbf{e}_j


最终得到无限均匀介质中集中脉冲力产生的静态位移场，

.. math:: 
    :label: Gtensor

    \mathbf{G} = \dfrac{1}{16 \pi \mu (1-\nu) r}  \left[ (3-4\nu)\mathbf{I} + \dfrac{\mathbf{rr^T}}{r^2} \right]

其中 :math:`\mathbf{G}` 为格林函数张量，其分量形式为

.. math:: 
    :label:

    G_{ij} = \dfrac{1}{16 \pi \mu (1-\nu) r}  \left[ (3-4\nu) \delta_{ij} + \gamma_i \gamma_j \right] = \dfrac{1}{4\pi \mu} \dfrac{\Delta}{1+\Delta} \dfrac{1}{R} \left( \dfrac{\delta_{ij}}{\Delta} + \gamma_i \gamma_j \right)

其中 :math:`\Delta=\dfrac{1}{3-4\nu}=\dfrac{\lambda+\mu}{\lambda+3\mu}` ， :math:`\gamma_i=r_i/r` 为方向余弦， :math:`r` 为源点和场点的直线距离。最终解 :eq:`Gtensor` 的表达式与朗道的弹性理论教材保持一致 [#]_ 。




------------------------------


.. [#] https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf
.. [#] https://en.wikipedia.org/wiki/Multidimensional_transform
.. [#] Л.Д. 朗道, Е.М. 栗弗席兹. 弹性理论[M]. 5版. 北京: 高等教育出版社, 2009:31.
