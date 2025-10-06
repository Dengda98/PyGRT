.. include:: common_OPTs.rst_


greenfn
=================

:简介: 使用广义反射透射系数矩阵法 + 离散波数法计算动态全波解格林函数

语法
-----------

**grt greenfn** 
|-M|\ *model*
|-D|\ *depsrc/deprcv*
|-N|\ *nt/dt*\ [**+w**\ *zeta*][**+n**\ *fac*]
|-R|\ *r1,r2,...*\ |\ *file*
|-O|\ *outdir*
[ |-H|\ *f1/f2* ]
[ |-L|\ *args* ]
[ |-K|\ [**+k**\ *k0*][**+s**\ *ampk*][**+e**\ *keps*][**+v**\ *vmin*] ]
[ |-E|\ *t0*\ [/*v0*] ]
[ |-P|\ *nthreads* ]
[ |-G|\ **e|v|h|s** ]
[ |-S|\ *i1,i2,...* ]
[ **-e** ]
[ **-s** ]
[ **-h** ]


描述
--------

**greenfn** 模块计算每个震中距 *r* 的格林函数波形保存路径为 
``{outdir}/{model}_{depsrc}_{deprcv}_{r}/{stype}.sac`` ，
其中支持的格林函数类型 ``stype`` 详见 :ref:`grn_types` 。
执行 **greenfn** 模块的原始命令会附加式地保存到 ``{outdir}/command`` 文件。

不同震源的格林函数单位为：

+ 爆炸源：  :math:`10^{-20} \, \frac{\text{cm}}{\text{dyne} \cdot \text{cm}}`  
+ 单力源：  :math:`10^{-15} \, \frac{\text{cm}}{\text{dyne}}`
+ 剪切源：  :math:`10^{-20} \, \frac{\text{cm}}{\text{dyne} \cdot \text{cm}}`



必选选项
----------

.. include:: explain_-M.rst_

.. include:: explain_-D.rst_

.. _-N:

**-N**\ *nt/dt*\ [**+w**\ *zeta*][**+n**\ *fac*]
    采样点数 *nt* 和采样时间间隔 *dt* (secs) ，这将决定计算的最高频。还可设置：

    + **+w**\ *zeta* - 虚频率系数 :math:`\zeta` [0.8]。Bouchon (1981) 提出在频率上添加一个虚频率偏移， 
      :math:`\omega \leftarrow \omega - i \zeta \dfrac{\pi}{T}` ，
      其中 :math:`T` 为时窗长度 nt\*dt，以使波数积分适当偏移实轴上的极点。
    
    + **+n**\ *fac* - 频率域插值倍数[1]。即在做逆傅里叶变换时，在频域上最高频后补零，
      相当于 *nt* ← *nt* \* *fac* , *dt* ← *dt* / *fac* ，使计算的波形更平滑。

    当时窗长度 nt\*dt 太小“包不住”有效信号，或时窗长度足够但时延(|-E|)不合适，输出的波形会发生混叠，
    此时需调整 |-N| 和 |-E| 。

.. include:: explain_-R.rst_

.. include:: explain_-O.rst_


可选选项
--------

.. _-H:

**-H**\ *f1/f2*
    待计算的频率范围（闭区间）(Hz) [默认计算由 |-D| 指定的全部范围]。 
    *f1* 和 *f2* 分别为最小最大频率，若取 -1 则分别对应低通和高通。

.. include:: explain_-L.rst_

.. _-K:

**-K**\ [**+k**\ *k0*][**+s**\ *ampk*][**+e**\ *keps*][**+v**\ *vmin*]
    控制波数积分上限

    .. math::

        k_{\text{max}} = \sqrt{ k_0 \cdot \dfrac{\pi}{\Delta h} + \textit{<ampk>} \cdot \left(\dfrac{\omega}{v_{\text{min}}}\right)^2}

    + **+k**\ *k0* - 控制零频的积分上限 [5.0]，其中深度差 :math:`\Delta h = \max(|z_s - z_r|, 1.0)` 。
    + **+s**\ *ampk* - 放大倍数 [1.15] 。
    + **+e**\ *keps* - 用于判断提前结束波数积分的收敛精度[0.0, 默认不使用]，详见
      Yao and Harkrider (1983) 和 :doc:`/Advanced/k_integ` 。
    + **+v**\ *vmin* - 参考最小速度，默认 :math:`\max{(\min\limits_{i} (\alpha_i \cup \beta_i), 0.1)}` 。  
      
      + 只要设置了 *vmin* ，不论正负， 是否启用快速收敛算法不再受 :math:`\Delta h` 控制 （见 |-D| ），
        相当于改为由 **+v**\ *vmin* 手动控制。
      + 当 *vmin* 为负数时，使用快速收敛算法（这个用法更多是用于 debug 测试）。

.. include:: explain_-E.rst_

.. include:: explain_-P.rst_

.. include:: explain_-G.rst_

.. _-S:

**-S**\ *i1,i2,...*
    指定频率索引值，输出对应频率下波数积分过程中的核函数文件，多个索引值用逗号分隔，文件保存目录为
    ``{outdir}_grtstats`` 。若使用 **-S**\ *-1* 则保存所有频点。
    关于文件格式及其读取详见 :doc:`/Advanced/integ_converg/integ_converg` 。

.. include:: explain_-egrn.rst_

.. include:: explain_-s.rst_

.. include:: explain_-h.rst_


示例
-------

详见教程：

+ :doc:`/Tutorial/dynamic/gfunc`
+ :doc:`/Advanced/integ_converg/integ_converg`
+ :doc:`/Advanced/kernel/kernel`