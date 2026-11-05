:author: 朱邓达
:date: 2026-11-01

.. include:: common_OPTs.rst_

synnoise
=================

:简介: 基于随机单力源激发的面波合成各台站的 ZNE 背景噪声速度记录

语法
-----------

**grt synnoise**
|-M|\ *model*
|-C|\ *sfile*
|-U|\ *rfile*
|-O|\ *dir*
[ |-N|\ [*n1*][/\ *n2*][/\ *dn*] ]
[ |-P|\ *nprocs* ]
[ **-s** ]
[ **-h** ]

描述
-----------

给定层状模型、震源位置和台站位置，自动计算 Rayleigh 和 Love 面波格林函数，
并合成 ZNE 三分量背景噪声速度记录，单位 **cm/s**。
每个震源按文件顺序等间隔激发一次，单力方向和强度在给定范围内随机抽取，
相同输入和随机种子可复现结果。

输出 SAC 文件的头段变量 **unused6**、 **unused7** 分别保存台站的北向、东向坐标，
单位 **km**。台站深度通过高程变量 **stel** 保存，其值为 **−1000 × depth**，单位 **m**。

必需选项
-----------

.. include:: explain_-M.rst_

.. _-C:

**-C**\ *sfile*
    震源文件。前五个有效行必须按以下顺序给出参数，随后每行为 **ID north east depth**：

    .. code-block:: text

       FORCE_Z        -1 1
       FORCE_AZIMUTH  0 360
       STRENGTH       1e15 2e15
       INTERVAL       0.5
       SEED           0

       # ID     north    east   depth
       SRC001       0      -1     0.0
       SRC002       1       2     0.5

    - **FORCE_Z**：单位力的向下分量范围，限于 **[−1, 1]**
    - **FORCE_AZIMUTH**：水平角范围，从北顺时针计，单位度，限于 **[0, 360]**
    - **STRENGTH**：非负单力强度范围，单位 **dyne**
    - **INTERVAL**：相邻震源的激发间隔，单位秒，必须为正；首个震源在 **0 s** 激发
    - **SEED**：随机种子，取非负整数

    三个范围均用空格分隔最小值和最大值，相等表示固定值。
    上例的力方向在整个球面上均匀抽取，源强度在给定范围内均匀抽取。
    坐标依次为北向、东向、向下深度，单位 **km**，深度非负。
    支持空行和 **#** 注释，ID 在文件内不可重复。
    激发时刻达到或超过记录长度的震源将跳过，超出记录末端的波形截断。
    不支持液体层内的单力源或零水平源台距。

.. _-U:

**-U**\ *rfile*
    台站文件。前两个有效行必须按以下顺序给出参数，随后每行为 **ID north east depth**：

    .. code-block:: text

       DT      0.1
       LENGTH  3600

       # ID     north    east   depth
       STA001       0       0     0.0
       STA002       0       1     0.5

    - **DT**：采样间隔，单位秒，必须为正
    - **LENGTH**：记录总长度，单位秒，必须为正，且 **LENGTH/DT** 必须为整数

    坐标约定与 |-C| 相同，支持空行和 **#** 注释。
    ID 用作输出台站目录名，在文件内不可重复，也不可包含路径分隔符。

.. _-O:

**-O**\ *dir*
    不存在或为空的输出目录。每个台站的完整记录保存为
    **dir/noise/台站ID/Z.sac**、**N.sac** 和 **E.sac**，起始时刻为 **0 s**。
    同时保存模型副本、频散文件 **phaseR.nc** 和 **phaseL.nc**、
    实际激发的震源信息 **sources.info**，以及台站文件副本 **receivers.info**。

输入、输出路径仅允许英文字母、数字及 **_ . - / :**，Windows 下还允许反斜杠。

可选选项
-----------

.. include:: explain_-Nmodes.rst_

.. _-P:

**-P**\ *nprocs*
    并行计算数，默认使用全部可用核心。

.. include:: explain_-silent.rst_
.. include:: explain_-h.rst_


示例
------------

+ :doc:`/Gallery/ex25/ex25`
