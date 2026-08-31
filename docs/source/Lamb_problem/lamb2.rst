:author: 朱邓达
:date: 2026-08-31

第二类 Lamb 问题
===================

第二类 Lamb 问题是指，在半空间模型中，震源位于地下、接收点位于地表的情形。
几何关系由震源到接收点的射线角和方位角描述，具体的参数、输出格式以及无量纲
结果恢复为物理量的关系，详见 :doc:`/Module/lamb2`。

.. tabs::

   .. group-tab:: CLI

      C 程序 :command:`grt` 提供了模块 :doc:`/Module/lamb2` 求解第二类 Lamb 问题。

      .. literalinclude:: run_lamb2/run.sh
         :language: bash
         :start-after: BEGIN LAMB2
         :end-before: END LAMB2

      使用重定向将结果保存到文件 *lamb2.txt* 中，其内容格式类似于

      .. literalinclude:: run_lamb2/head_lamb2
         :language: text

   .. group-tab:: Python

      Python 提供了函数 :func:`lamb2() <pygrt.utils.lamb2>`，一次返回 Green 函数、
      相对于源点坐标的一阶空间偏导和相对于接收点坐标的一阶空间偏导。

      .. literalinclude:: run_lamb2/lamb2_plot_time.py
         :language: python
         :start-after: BEGIN LAMB2
         :end-before: END LAMB2


最后绘制计算得到的格林函数以及相对源点坐标的一阶空间导数。

:download:`lamb2_plot_time.py <run_lamb2/lamb2_plot_time.py>`

.. figure:: run_lamb2/lamb2.svg
   :align: center

   复现了原书中的图 7.4.6

-----------

.. figure:: run_lamb2/lamb2_d1.svg
   :align: center

   复现了原书中的图 7.4.10

-----------

.. figure:: run_lamb2/lamb2_d2.svg
   :align: center

   复现了原书中的图 7.4.11

-----------

.. figure:: run_lamb2/lamb2_d3.svg
   :align: center

   复现了原书中的图 7.4.12

频域解和时域解的对比
-------------------------------

由于 PyGRT 的频域解中可以计算格林函数相对于接收点坐标的空间偏导，
因此这里我们可以进行更多的对比。在以下对比图中发现，对于频域解卷积了阶跃函数之后，
格林函数的 Gibbs 效应少了很多，然而对于格林函数的空间导数还是很明显，
这是因为经过理论推导，空间导数项转为了时间导数项。
如果想要更清晰的对比，可以自行对空间导数项再做一次积分然后绘制。

:download:`lamb2_plot_freq_time.py <run_lamb2/lamb2_plot_freq_time.py>`

.. figure:: run_lamb2/lamb2_compare_freq_time.svg
   :align: center

-----------

.. figure:: run_lamb2/lamb2_compare_freq_time_z.svg
   :align: center

-----------

.. figure:: run_lamb2/lamb2_compare_freq_time_r.svg
   :align: center


关于格林函数对接收点坐标的一阶空间偏导的数学推导
---------------------------------------------------

记号和归一化
~~~~~~~~~~~~~~

设接收点坐标为 :math:`\boldsymbol{x}=(x_1,x_2,x_3)`，源点坐标为
:math:`\boldsymbol{x}'=(x'_1,x'_2,x'_3)`。下标 :math:`i` 表示位移分量，
:math:`j` 表示点力方向，:math:`a` 表示接收点坐标方向，:math:`a'` 表示源点
坐标方向，其中 :math:`a=1,2,3`。为避免与源点坐标方向的撇号混淆，波速比及其
互补量记为

.. math::

   \kappa=\frac{\beta}{\alpha},\qquad
   \kappa'=\sqrt{1-\kappa^2},\qquad
   \chi=1-2\kappa^2

这里的 :math:`\kappa'` 是材料参数的互补量，不表示源点坐标导数。空间偏导的
逗号记号定义为

.. math::

   \bar G^H_{ij,a}
   \equiv
   \frac{\partial \bar G^H_{ij}}{\partial x_a},
   \qquad
   \bar G^H_{ij,a'}
   \equiv
   \frac{\partial \bar G^H_{ij}}{\partial x'_a}.

第二类 Lamb 问题的阶跃力 Green 函数及其两类空间偏导，均由 P 项、S 项和
S-P 项的时间函数组合得到。按照第 7.1.2 节以及第 7.4.1.4 节的组合关系，
接收点导数为

.. math::

   \bar G^H_{ij,a}(\bar t)
   =
   \frac{\partial}{\partial\bar t}
   \left[
   F_{P,ij,a}(\bar t)
   +
   F_{S,ij,a}(\bar t)
   -
   H(\theta-\theta_c)F_{S-P,ij,a}(\bar t)
   \right],
   \qquad
   \theta_c=\arcsin\kappa

源点导数具有完全相同的时间组合形式

.. math::

   \bar G^H_{ij,a'}(\bar t)
   =
   \frac{\partial}{\partial\bar t}
   \left[
   F_{P,ij,a'}(\bar t)
   +
   F_{S,ij,a'}(\bar t)
   -
   H(\theta-\theta_c)F_{S-P,ij,a'}(\bar t)
   \right].

其中 :math:`F_{S-P,ij,a}` 和 :math:`F_{S-P,ij,a'}` 仅在
:math:`\bar t_{S-P}<\bar t<1` 的相应时间区间内参与组合。若 :math:`G^H` 表示
具有物理量纲的阶跃力 Green 函数，则无量纲量与物理量之间的关系为

.. math::

   \bar G^H_{ij}
   =
   \pi^2\mu rG^H_{ij},
   \qquad
   \bar G^H_{ij,a'}
   =
   \pi^2\mu r^2G^H_{ij,a'},
   \qquad
   \bar G^H_{ij,a}
   =
   \pi^2\mu r^2G^H_{ij,a}.

因此 Green 函数的恢复因子为 :math:`\pi^2\mu r`，两类空间偏导的恢复因子
均为 :math:`\pi^2\mu r^2`

源点与接收点导数的关系
~~~~~~~~~~~~~~~~~~~~~~~~~~

水平坐标方向
^^^^^^^^^^^^

均匀半空间对水平平移不变，积分核中的水平坐标只以
:math:`x_a-x'_a` 的形式出现。因此，对 :math:`a=1,2` 有

.. math::

   \frac{\partial}{\partial x_a}
   =
   -\frac{\partial}{\partial x'_a},
   \qquad
   \bar G^H_{ij,a}=-\bar G^H_{ij,a'},
   \qquad a=1,2.

同一关系对 P、S 和 S-P 三项的矩阵分子分别成立。若 :math:`M` 表示 P 项
矩阵，:math:`N` 表示 S 项及 S-P 项矩阵，则

.. math::

   M_{ij,a}=-M_{ij,a'},
   \qquad
   N_{ij,a}=-N_{ij,a'},
   \qquad a=1,2.

接收点竖直方向
^^^^^^^^^^^^^^^^

接收点位于自由表面时，位移 Green 函数列
:math:`u_i=G^H_{ij}` 满足自由表面牵引条件

.. math::

   \sigma_{13}
   =\mu\left(u_{1,3}+u_{3,1}\right)=0,
   \qquad
   \sigma_{23}
   =\mu\left(u_{2,3}+u_{3,2}\right)=0,

.. math::

   \sigma_{33}
   =
   \lambda\left(u_{1,1}+u_{2,2}\right)
   +(\lambda+2\mu)u_{3,3}
   =0.

于是

.. math::

   u_{1,3}=-u_{3,1},
   \qquad
   u_{2,3}=-u_{3,2},
   \qquad
   u_{3,3}
   =
   -\frac{\lambda}{\lambda+2\mu}
   \left(u_{1,1}+u_{2,2}\right),
   \qquad
   \frac{\lambda}{\lambda+2\mu}=\chi.

结合水平坐标方向的平移关系，接收点竖直导数的前两行可由源点水平导数
给出。对 Green 函数矩阵可写为

.. math::

   \left[\bar G^H_{ij,3}\right]_{i,j=1}^{3}
   =
   \begin{bmatrix}
   \bar G^H_{31,1'} & \bar G^H_{32,1'} & \bar G^H_{33,1'} \\
   \bar G^H_{31,2'} & \bar G^H_{32,2'} & \bar G^H_{33,2'} \\
   \bar G^H_{31,3}  & \bar G^H_{32,3}  & \bar G^H_{33,3}
   \end{bmatrix}.

对每一类积分矩阵，同样有

.. math::

   M_{,3}
   =
   \begin{bmatrix}
   M_{31,1'} & M_{32,1'} & M_{33,1'} \\
   M_{31,2'} & M_{32,2'} & M_{33,2'} \\
   M_{31,3}  & M_{32,3}  & M_{33,3}
   \end{bmatrix},
   \qquad
   N_{,3}
   =
   \begin{bmatrix}
   N_{31,1'} & N_{32,1'} & N_{33,1'} \\
   N_{31,2'} & N_{32,2'} & N_{33,2'} \\
   N_{31,3}  & N_{32,3}  & N_{33,3}
   \end{bmatrix}.

利用源点导数矩阵的对称关系，还可以使用
:math:`M_{31,2'}=M_{32,1'}` 和 :math:`N_{31,2'}=N_{32,1'}` 的形式。上述矩阵
的最后一行是新的接收点竖直导数，需要由自由表面关系和相应的积分分母
有理化得到

P 波项的分子
~~~~~~~~~~~~~~~~

按照第 7.2 节，将 P 波项的变量取为 :math:`x=\eta_\alpha`，并定义

.. math::

   g=x^2-\kappa^2,
   \qquad
   \eta_\beta^2=x^2+\kappa'^2,
   \qquad
   \gamma_P=2x^2+\chi,
   \qquad
   q=\frac{x\cos\theta-\bar t}{\sin\theta},

.. math::

   p^2
   =
   \frac{x^2-2x\bar t\cos\theta+\bar t^2-\kappa^2\sin^2\theta}
   {\sin^2\theta}.

式 (7.2.3a)--(7.2.3b) 给出的分母关系为

.. math::

   \frac{1}{R^P(p,q)}
   =
   \frac{\gamma_P^2+4x\eta_\beta g}{R'^P(x)}.

在此分母约定下，各项分子写成

.. math::

   \frac{M^{P(1)}_{ij,a}(x)}{R'^P(x)}
   +
   \frac{M^{P(2)}_{ij,a}(x)}
   {R'^P(x)\eta_\beta}.

源点竖直方向的一阶导数由式 (7.1.6c) 给出。在上述变量变换下，对
:math:`ξ=1,2` 有

.. math::

   M^{P(\xi)}_{ij,3'}(x)
   =
   -xM^{P(\xi)}_{ij}(x).

在尚未进行分母拆分时，P 波矩阵最后一行为

.. math::

   \begin{aligned}
   M^P_{31,3}
   &=-2\chi q\eta_\beta(q^2-p^2)\cos\phi, &
   M^P_{32,3}
   &=-2\chi q\eta_\beta(q^2-p^2)\sin\phi, \\
   M^P_{33,3}
   &=-2\chi\eta_\alpha\eta_\beta(q^2-p^2).
   \end{aligned}

令 :math:`\xi=1,2`，由接收点竖直导数前两行与源点水平导数的关系，有

.. math::

   \begin{aligned}
   M^{P(\xi)}_{11,3}&=M^{P(\xi)}_{31,1'},&
   M^{P(\xi)}_{12,3}&=M^{P(\xi)}_{32,1'},&
   M^{P(\xi)}_{13,3}&=M^{P(\xi)}_{33,1'},\\
   M^{P(\xi)}_{21,3}&=M^{P(\xi)}_{31,2'},&
   M^{P(\xi)}_{22,3}&=M^{P(\xi)}_{32,2'},&
   M^{P(\xi)}_{23,3}&=M^{P(\xi)}_{33,2'}.
   \end{aligned}

将上面的最后一行代入式 (7.2.3a)--(7.2.3b) 并按 :math:`R'^P(x)` 和
:math:`R'^P(x)\eta_\beta` 拆分，得到

.. math::

   \begin{aligned}
   M^{P(1)}_{31,3}
   &=-8\chi x^2g\eta_\beta^2q(q^2-p^2)\cos\phi,\\
   M^{P(1)}_{32,3}
   &=-8\chi x^2g\eta_\beta^2q(q^2-p^2)\sin\phi,\\
   M^{P(1)}_{33,3}
   &=-8\chi x^3g\eta_\beta^2(q^2-p^2),\\[1mm]
   M^{P(2)}_{31,3}
   &=-2\chi xq\eta_\beta^2\gamma_P^2(q^2-p^2)\cos\phi,\\
   M^{P(2)}_{32,3}
   &=-2\chi xq\eta_\beta^2\gamma_P^2(q^2-p^2)\sin\phi,\\
   M^{P(2)}_{33,3}
   &=-2\chi x^2\eta_\beta^2\gamma_P^2(q^2-p^2).
   \end{aligned}

S 波项和 S-P 项的分子
~~~~~~~~~~~~~~~~~~~~~~~~~~

按照第 7.3 节，将 S 波项的变量取为 :math:`x=\eta_\beta`，并定义

.. math::

   g'=x^2-1,
   \qquad
   \eta_\alpha^2=x^2-\kappa'^2,
   \qquad
   \gamma_S=2x^2-1,
   \qquad
   q=\frac{x\cos\theta-\bar t}{\sin\theta},

.. math::

   p^2
   =
   \frac{x^2-2x\bar t\cos\theta+\bar t^2-\sin^2\theta}
   {\sin^2\theta}.

式 (7.3.1a)--(7.3.1f) 中的分母按第 7.3 节有理化为

.. math::

   \frac{1}{R^S(p,q)}
   =
   \frac{\gamma_S^2+4x\eta_\alpha g'}{R'^S(x)}.

相应的分子写成

.. math::

   \frac{N^{S(1)}_{ij,a}(x)}{R'^S(x)}
   +
   \frac{N^{S(2)}_{ij,a}(x)}
   {R'^S(x)\eta_\alpha}.

源点竖直方向的一阶导数在第 7.3 节的变量变换下满足

.. math::

   N^{S(\xi)}_{ij,3'}(x)
   =
   -xN^{S(\xi)}_{ij}(x),
   \qquad \xi=1,2.

在尚未进行分母拆分时，S 波矩阵最后一行为

.. math::

   \begin{aligned}
   N^S_{31,3}
   &=-\chi q\eta_\beta\gamma_S\cos\phi, &
   N^S_{32,3}
   &=-\chi q\eta_\beta\gamma_S\sin\phi, \\
   N^S_{33,3}
   &=\chi\gamma_S(q^2-p^2).
   \end{aligned}

对 :math:`\xi=1,2`，其前两行满足

.. math::

   \begin{aligned}
   N^{S(\xi)}_{11,3}&=N^{S(\xi)}_{31,1'},&
   N^{S(\xi)}_{12,3}&=N^{S(\xi)}_{32,1'},&
   N^{S(\xi)}_{13,3}&=N^{S(\xi)}_{33,1'},\\
   N^{S(\xi)}_{21,3}&=N^{S(\xi)}_{31,2'},&
   N^{S(\xi)}_{22,3}&=N^{S(\xi)}_{32,2'},&
   N^{S(\xi)}_{23,3}&=N^{S(\xi)}_{33,2'}.
   \end{aligned}

将最后一行代入式 (7.3.1a)--(7.3.1f) 并完成分母拆分，可得

.. math::

   \begin{aligned}
   N^{S(1)}_{31,3}
   &=-\chi x^2q\gamma_S^3\cos\phi,\\
   N^{S(1)}_{32,3}
   &=-\chi x^2q\gamma_S^3\sin\phi,\\
   N^{S(1)}_{33,3}
   &=\chi x\gamma_S^3(q^2-p^2),\\[1mm]
   N^{S(2)}_{31,3}
   &=-4\chi x^3g'\eta_\alpha^2\gamma_Sq\cos\phi,\\
   N^{S(2)}_{32,3}
   &=-4\chi x^3g'\eta_\alpha^2\gamma_Sq\sin\phi,\\
   N^{S(2)}_{33,3}
   &=4\chi x^2g'\eta_\alpha^2\gamma_S(q^2-p^2).
   \end{aligned}

上述两组 :math:`N^{S(\xi)}_{ij,a}` 同时用于 S1、S2 和 S-P 三类积分，
三者的差异只来自第 7.3.5 节和式 (7.3.9) 所规定的积分路径、积分上下限
以及基本积分的取值。时间函数的组合关系为

.. math::

   F_{S,ij,a}
   =
   F_{S1,ij,a}
   -
   H(\theta-\theta_c)F_{S2,ij,a},
   \qquad
   F_{S,ij,a'}
   =
   F_{S1,ij,a'}
   -
   H(\theta-\theta_c)F_{S2,ij,a'}.

公式索引
~~~~~~~~

上述关系所依据的书中公式可按求解顺序归纳如下

* 基本 P、S 矩阵：第 7.1.2 节式 (7.1.3a)--(7.1.3b)
* 源点坐标的一阶空间导数：式 (7.1.4)--(7.1.6f)
* P 波项的变量、矩阵分子和分母：式 (7.2.1)--(7.2.3b)
* S 波项及 S-P 项的变量、矩阵分子和分母：式 (7.3.1a)--(7.3.1f)
* S1、S2 与 S-P 项的积分路径和组合：第 7.3.5 节及式 (7.3.9)
* 广义闭合解的时间组合与导数：第 7.4 节，尤其是第 7.4.1.4 节
