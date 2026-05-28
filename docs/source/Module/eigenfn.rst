:author: 朱邓达
:date: 2026-05-20

.. include:: common_OPTs.rst_


eigenfn
=================

:简介: 基于面波频散结果，计算面波本征函数、能量积分、群速度、敏感核

语法
-----------

**grt eigenfn** 
|-C|\ *path*
[ |-F|\ *f1*\ [/\ *f2*][/\ *df*][**+p**] ]
[ |-N|\ [*n1*][/\ *n2*][/\ *dn*] ]
[ |-W|\ *path*\ [**+z**\ *z1*\ [/\ *z2/dz*]] ]
[ |-U|\ *path* ]
[ |-K|\ [**+c**\ *path*][**+u**\ *path*][**+x**\ *path*][**+z**\ *dz*] ]
[ **-h** ]

描述
--------

基于 **eigenv** 模块计算的相速度频散结果， **以及其结果文件中保存的模型路径，** 
**eigenfn** 模块可计算出多个物理量，各个量之间的依赖关系为

.. container:: small-mermaid

    .. mermaid::

        flowchart LR 

        V(["Phase-Velocity Dispersion<br/> (Eigenvalues)"])
        F(["Eigenfunctions"])
        E(["Energy Integrals"])
        G(["Group-Velocity Dispersion"])
        K(["Phase/Group-Velocity Sensitivity"])

        V --> F 
        F --> E 
        E --> G 
        E --> K

        classDef cmdcls fill:#f9f2d9,stroke:#e8d174,stroke-width:2px,color:#333;
        class V cmdcls

各个量的结果均以 |NetCDF| 网格格式保存，涉及的选项如下，其中公式符号的定义详见 |yao2026p| ：

+ **本征函数（Eigenfunctions）**： |-W|
  
  ``ncdump -h`` 命令查看本征函数结果，输出类似于：

  ::

      netcdf egn_R {
      dimensions:
              freq = 1 ;
              mode = 11 ;
              z = 61 ;
              w = 8 ;
      variables:
              double freq(freq) ;
              int mode(mode) ;
              double z(z) ;
              int cnum(freq) ;
              double c(freq, mode) ;
                      c:_FillValue = 0. ;
              double eigfn(freq, mode, z, w) ;
                      eigfn:_FillValue = 0. ;

      // global attributes:
                      :isRayl = 1 ;
      }

  .. grid:: 3
      :gutter: 2

      .. grid-item:: **维度（dimensions）**
          :columns: auto 

          + ``freq`` - 频率点数
          + ``mode`` - 阶数点数
          + ``z`` - 深度点数
          + ``w`` - 本征函数个数 [#]_
      
      .. grid-item:: **变量（variables）**
          :columns: auto 

          + ``freq`` - 频率数组 (Hz)
          + ``mode`` - 阶数数组
          + ``z`` - 深度数组 (km)
          + ``cnum`` - 不同频率下的相速度点数
          + ``c`` - 不同频率不同阶的相速度 (km/s)
          + ``eigfn`` - 不同频率不同阶不同深度的本征函数
      
      .. grid-item:: **全局属性（global attributes）**
          :columns: auto 

          + ``isRayl`` - 是否为Rayleigh波结果

  .. [#] ``w = 8`` 对应 Rayleigh 波，分别为
    :math:`\text{Re}[q],\text{Im}[q],\text{Re}[w],\text{Im}[w],\text{Re}[\sigma_R],\text{Im}[\sigma_R],\text{Re}[\tau_R],\text{Im}[\tau_R]` ，
    ``w = 4`` 对应 Love 波，分别为
    :math:`\text{Re}[v],\text{Im}[v],\text{Re}[\tau_L],\text{Im}[\tau_L]`    

+ **群速度频散（Group-Velocity Dispersion）**： |-U|

  ``ncdump -h`` 命令查看群速度频散结果，输出类似于：

  ::

      netcdf group_R {
      dimensions:
              freq = 1000 ;
              mode = 49 ;
      variables:
              double freq(freq) ;
              int mode(mode) ;
              int cnum(freq) ;
              double c(freq, mode) ;
                      c:_FillValue = 0. ;
              double u(freq, mode) ;
                      u:_FillValue = 0. ;

      // global attributes:
                      :isRayl = 1 ;
      }

  .. grid:: 3
      :gutter: 2

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
          + ``u`` - 不同频率不同阶的群速度 (km/s)
      
      .. grid-item:: **全局属性（global attributes）**
          :columns: auto 

          + ``isRayl`` - 是否为Rayleigh波结果


+ **能量积分（Energy Integrals）**：|-K|\ **+x**
  
  ``ncdump -h`` 命令查看能量积分结果，输出类似于：

  ::

      netcdf egy_R {
      dimensions:
              freq = 500 ;
              mode = 1 ;
              e = 10 ;
      variables:
              double freq(freq) ;
              int mode(mode) ;
              int cnum(freq) ;
              double c(freq, mode) ;
                      c:_FillValue = 0. ;
              double egyint(freq, mode, e) ;
                      egyint:_FillValue = 0. ;

      // global attributes:
                      :isRayl = 1 ;
      }

  .. grid:: 3
      :gutter: 2

      .. grid-item:: **维度（dimensions）**
          :columns: auto 

          + ``freq`` - 频率点数
          + ``mode`` - 阶数点数
          + ``e`` - 能量积分个数 [#]_
      
      .. grid-item:: **变量（variables）**
          :columns: auto 

          + ``freq`` - 频率数组 (Hz)
          + ``mode`` - 阶数数组
          + ``cnum`` - 不同频率下的相速度点数
          + ``c`` - 不同频率不同阶的相速度 (km/s)
          + ``egyint`` - 能量积分结果
      
      .. grid-item:: **全局属性（global attributes）**
          :columns: auto 

          + ``isRayl`` - 是否为Rayleigh波结果

  .. [#] 能量积分结果始终为 10 个，分别对应 
      :math:`\text{Re}[I_1],\text{Im}[I_1],\text{Re}[I_2],\text{Im}[I_2],\text{Re}[I_3],\text{Im}[I_3],\text{Re}[I_4],\text{Im}[I_4],\text{Re}[\varepsilon],\text{Im}[\varepsilon]`，
      其中 :math:`\varepsilon = \omega^2 I_1 - k^2 I_2 + k I_3 - I_4` 为能量积分的验证项，应极接近 0。对于 Love 波， :math:`I_3` 为 0 。

+ **相/群速度敏感核（Phase/Group-velocity Sensitivity）**: |-K|\ **+c** , |-K|\ **+u**

  ``ncdump -h`` 命令查看相/群速度敏感核结果，输出类似于：

  ::

      netcdf csens_R {
      dimensions:
              freq = 500 ;
              mode = 1 ;
              z = 31 ;
              k = 3 ;
      variables:
              double freq(freq) ;
              int mode(mode) ;
              double z(z) ;
              int cnum(freq) ;
              double c(freq, mode) ;
                      c:_FillValue = 0. ;
              double csens(freq, mode, z, k) ;
                      csens:_FillValue = 0. ;

      // global attributes:
                      :isRayl = 1 ;
      }

  .. grid:: 3
      :gutter: 2

      .. grid-item:: **维度（dimensions）**
          :columns: auto 

          + ``freq`` - 频率点数
          + ``mode`` - 阶数点数
          + ``z`` - 深度点数
          + ``k`` - 敏感核个数 [#]_
      
      .. grid-item:: **变量（variables）**
          :columns: auto 

          + ``freq`` - 频率数组 (Hz)
          + ``mode`` - 阶数数组
          + ``z`` - 深度数组 (km)
          + ``cnum`` - 不同频率下的相速度点数

          对于相速度敏感核，变量包括：

          + ``c`` - 不同频率不同阶的相速度 (km/s)
          + ``csens`` - 不同频率不同阶不同深度的相速度敏感核
          
          对于群速度敏感核，变量包括：

          + ``u`` - 不同频率不同阶的群速度 (km/s)
          + ``usens`` - 不同频率不同阶不同深度的群速度敏感核
      
      .. grid-item:: **全局属性（global attributes）**
          :columns: auto 

          + ``isRayl`` - 是否为Rayleigh波结果

  .. [#] 敏感核个数始终为3个，分别对应相/群速度对于 :math:`\alpha,\beta,\rho` 的无量纲化偏导；
      对于 Love 波结果，相/群速度对于 :math:`\alpha` 的偏导为 0 。


必选选项
----------

.. include:: explain_-Ceigv.rst_


可选选项
--------

.. include:: explain_-Ffreqs.rst_

.. include:: explain_-Nmodes.rst_

.. _-W:

**-W**\ *path*\ [**+z**\ *z1*\ [/\ *z2/dz*]] 
    计算不同深度上的本征函数，以 |NetCDF| 格式输出到路径 *path* 。
    关于 **+z**，有以下几种用法：

    + 如果不指定 **+z**，则计算模型界面深度上的本征函数；
    + **+z**\ *z1* - 单点深度 *z1* ，单位 km
    + **+z**\ *z1/z2/dz* - 从深度 *z1* 到 *z2* ，间距 *dz* 的一系列等距深度，单位 km
    
.. _-U:

**-U**\ *path*
    计算相速度频散对应的群速度频散，以 |NetCDF| 格式输出到路径 *path* 。

.. _-K:

**-K**\ [**+c**\ *path*][**+u**\ *path*][**+x**\ *path*][**+z**\ *dz*]
    计算能量积分和相/群速度敏感核，以 |NetCDF| 格式输出到路径 *path* 。

    + **+c**\ *path* - 输出相速度敏感核（无量纲化）
      :math:`\dfrac{\alpha}{c}\dfrac{\partial c}{\partial \alpha}, \dfrac{\beta}{c}\dfrac{\partial c}{\partial \beta}, \dfrac{\rho}{c}\dfrac{\partial c}{\partial \rho}`
    + **+u**\ *path* - 输出群速度敏感核（无量纲化）
      :math:`\dfrac{\alpha}{U}\dfrac{\partial U}{\partial \alpha}, \dfrac{\beta}{U}\dfrac{\partial U}{\partial \beta}, \dfrac{\rho}{U}\dfrac{\partial U}{\partial \rho}`
    + **+x**\ *path* - 输出能量积分 :math:`I_1,I_2,I_3,I_4` 
      以及验证项 :math:`\omega^2 I_1 - k^2 I_2 + k I_3 - I_4` 
    + **+z**\ *dz* - 计算能量积分和相/群速度敏感核时使用的深度间隔（km），
      这相当于将层状模型进行阶梯式插值后再进行计算 [默认直接基于模型各层厚度]。


.. include:: explain_-h.rst_


示例
-------
