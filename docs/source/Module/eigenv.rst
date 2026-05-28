:author: 朱邓达
:date: 2026-05-20

.. include:: common_OPTs.rst_


eigenv
=================

:简介: 计算面波频散曲线

语法
-----------

**grt eigenv** 
|-M|\ *model*
|-C|\ *path*
|-S|\ **R|L**
|-F|\ *f1*\ [/\ *f2/df*][**+p**]
[ |-N|\ [*max_order*] ]
[ |-X|\ *freq*\ [**+c**\ *cmin*/\ *cmax*][**+i**\ *iref*] ]
[ |-T|\ [**+t**\ *tol*][**+c**\ *cgap*][**+r**\ *thrd*][**+v**\ *vgap*][**+u**\ *dc*] ]
[ |-P|\ *nthreads* ]
[ **-s** ]
[ **-h** ]


描述
--------

通过搜索由广义反射透射系数矩阵法定义的久期函数的零点来计算层状介质的面波相速度频散曲线，
结果保存为 |NetCDF| 网格格式。
如果使用 ``ncdump -h`` 命令查看该文件的基本信息，结果类似于：
::

    netcdf phase_R {
    dimensions:
            freq = 500 ;
            mode = 54 ;
    variables:
            double freq(freq) ;
            int mode(mode) ;
            int cnum(freq) ;
            double c(freq, mode) ;
                    c:_FillValue = 0. ;
            int ciref(freq, mode) ;
                    ciref:_FillValue = -1 ;

    // global attributes:
                    :command = "grt eigenv -Mmod2 -SR -F0/5/0.01 -N -Cphase_R.nc" ;
                    :model = "mod2" ;
                    :isRayl = 1 ;
    }

.. grid:: 3
    :gutter: 4
    :margin: 3

    .. grid-item:: **维度（dimensions）**
        :columns: auto 

        + ``freq`` - 频率点数
        + ``mode`` - 阶数点数
    
    .. grid-item:: **变量（variables）**
        :columns: auto 

        + ``freq`` - 频率数组 (Hz)
        + ``mode`` - 阶数数组
        + ``cnum`` - 不同频率下的相速度点数
        + ``c`` - 不同频率不同阶的相速度 (km/s)
        + ``ciref`` - 搜索对应相速度使用的久期函数层位
    
    .. grid-item:: **全局属性（global attributes）**
        :columns: auto 

        + ``command`` - 生成该文件的命令
        + ``model`` - 模型路径
        + ``isRayl`` - 是否为Rayleigh波结果


必选选项
----------

.. include:: explain_-M.rst_

.. include:: explain_-Ceigv.rst_

.. _-S:

**-S**\ **R|L**
    计算 **R**\ ayleigh 波频散还是 **L**\ ove 波频散。

.. _-F:

**-F**\ *f1*\ [/\ *f2/df*][**+p**]
    设置要计算的频率，有以下几种用法：

    + **-F**\ *f1* - 单个频率 *f1* ，单位 Hz
    + **-F**\ *f1/f2/df* - 从频率 *f1* 到 *f2* ，间距 *df* 的一系列等距频率，单位 Hz
    
    在以上用法后加上 **+p** ，则 *f1,f2,df* 均变为周期，单位 secs。
    例如 **-F**\ 0/10/2+p 表示计算周期 2s, 4s, 6s, 8s, 10s 处的频散。
    
    **注：** 在使用 *f1/f2/df* 定义等距频率点时，为了方便，用户可以设置 *f1=0*，
    此时程序会自动跳到下一个频率点 (*df*)。


可选选项
--------

.. _-N:

**-N**\ [*max_order*]
    设置计算的最大阶数，0 对应基阶。

    + 如果不使用 **-N** ，默认仅计算基阶，相当于设置 **-N**\ *0*
    + 如果省略 *max_order* ，即仅给定 **-N** ，则尽可能计算所有阶。

.. _-X:

**-X**\ *freq*\ [**+c**\ *cmin*/\ *cmax*][**+i**\ *iref*]
    该选项常用于 debug，可输出频率为 *freq* Hz，相速度范围在 [\ *cmin*, *cmax*] 之间的第 *iref* 层的久期函数。
    此时 **-F** 和 **-C** 不必再设置，久期函数结果将输出到标准输出，包含三列，分别为相速度、久期函数实部、久期函数虚部。
    采样间隔受搜根方法控制，详见 |-T| 。

.. _-T:

**-T**\ [**+t**\ *tol*][**+c**\ *cgap*][**+r**\ *thrd*][**+v**\ *vgap*][**+u**\ *dc*]
    设置搜根过程中的一些参数，绝大多数情况下默认参数即可。
    默认情况下， **eigenv** 模块使用自适应采样久期函数的方式进行搜根 (|ars2026|)，
    其依赖参数：

    + **+t**\ *tol* - 自适应策略的收敛容差 :math:`10^{-\it \text{tol}}` ， *tol* 默认为 3。
    
    作为测试，也可使用等相速度间隔的方式进行搜根，其依赖参数：

    + **+u**\ *dc* - 等间隔搜索的相速度间隔，单位 km/s 。

    搜根过程也依赖以下参数：

    + **+c**\ *cgap* - 同频率下相邻零点的最小相速度间隔 :math:`10^{-\it \text{cgap}}` ， *cgap* 默认为 7。
    + **+r**\ *thrd* - 识别为零点的久期函数振幅阈值 :math:`10^{-\it \text{thrd}}` ， *thrd* 默认为 4。
    + **+v**\ *vgap* - 零点位置与模型速度的最小间隔 :math:`10^{-\it \text{vgap}}` ， *vgap* 默认为 7。
    
.. include:: explain_-P.rst_
.. 
    
    **eigenv** 模块会自动并行搜索不同频率的频散结果，提高计算效率。

.. include:: explain_-s.rst_

.. include:: explain_-h.rst_


示例
-------


