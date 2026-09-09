:author: 朱邓达
:date: 2026-09-08

.. include:: common_OPTs.rst_


rcvfn
=================

:简介: 计算平面入射 P 或 SV 波的接收函数及自由表面位移响应

语法
-----------

**grt rcvfn**
|-M|\ *model*
( |-P|\ *rayp* | |-I|\ *inca*\ [/\ *idx*] )
|-T|\ **P|S**
|-N|\ *nt/dt*\ [**+w**\ *zeta*][**+n**\ *fac*][**+a**][**+f**]
|-O|\ *outdir*
[ |-A|\ *alp* ]
[ |-E|\ *delay* ]
[ |-W| ]
[ **-s** ]
[ **-h** ]


描述
--------

**rcvfn** 将入射波视为沿指定水平射线参数传播的单位平面波，在一维水平层状半空间中计算其到达自由表面后的 P-SV 响应。
入射波从模型底部半空间，或 |-I| 指定的模型层，向上入射。

默认只保存接收函数；设置 |-W| 后，还会保存单位入射位移对应的 Z/R 响应。
具体公式详见 :doc:`/Formula/others/rcvfn`。

SAC 头段中的 ``user0``、``user1``、``user2`` 分别记录虚频系数 :math:`\omega_I`、水平射线参数 *rayp* 和高斯滤波参数 *alp*。


必选选项
----------

.. include:: explain_-M.rst_

.. _-P:

**-P**\ *rayp*
    直接设置水平射线参数，单位为 s/km，要求为正且有限。它与 |-I| 互斥。
    对于入射波速度 :math:`V` 和从竖直方向量起的入射角 :math:`\theta`，二者满足
    :math:`\mathit{rayp}=\sin\theta/V`。实际入射波类型由 **-T** 指定。

.. _-I:

**-I**\ *inca*\ [/\ *idx*]
    通过入射角设置水平射线参数，入射角单位为度，取值范围为 :math:`[0,90)`。
    *idx* 是从模型底部反向计数的层号：*idx=0* 表示底部半空间，``idx=1`` 表示其上方一层，依此类推；省略时默认为 ``idx=0``。
    程序使用所选层的入射波速度按 :math:`\mathit{rayp}=\sin\theta/V` 换算，并将模型截取到该入射层。
    它与 **-P** 互斥。

.. _-T:

**-T**\ **P|S**
    选择入射波类型。``P`` 表示 P 波，``S`` 表示 SV 波；本模块的 ``S`` 不表示 SH 波。

.. include:: explain_-Nnt.rst_

.. include:: explain_-O.rst_


可选选项
--------

.. _-A:

**-A**\ *alp*
    设置高斯低通滤波参数 *alp*，单位为 Hz，默认值为 1.0。频域滤波函数为

    .. math::

        H(f)=\exp\left[-\left(\frac{\pi f}{\mathit{alp}}\right)^2\right]

    对应的特征拐角频率约为 :math:`\mathit{alp}/\pi`。

.. _-E:

**-E**\ *delay*
    增加附加延迟，单位为 s。 **程序内部会自动安排一部分延迟**，
    只有在需要更长的起始时间时才需要设置该选项。
    尤其是 S 波入射时，S 波到时前后都需要足够的时窗以防止混叠，此时需调整 |-N| 和 |-E| 。

.. _-W:

**-W**
    除接收函数外，同时输出单位入射位移对应的垂向和径向响应。

.. include:: explain_-silent.rst_

.. include:: explain_-h.rst_

示例
-------

+ :doc:`/Tutorial/receiver_function/rcvfn`
