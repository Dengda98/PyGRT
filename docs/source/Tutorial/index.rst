简易教程
=============

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

**PyGRT** 程序包由C和Python两个编程语言的代码组成，目的是兼并高效性和便捷性。底层复杂运算由C语言编写，编译链接成动态库 ``libgrt.so`` 供Python调用。Python通过 ``ctypes`` 库导入动态库以使用外部函数，以此兼并了C语言的高效和Python语言的便捷。

除了Python脚本式运行， **PyGRT** 保留传统命令行式运行C程序，即在建立动态库过程中也编译了一些可执行文件。 **C程序的运行独立于Python，不需要Python环境，以满足更多计算场景。** 每个C程序可使用 ``-h`` 查看帮助。

计算 **动态解** 的主要计算流程如下：

+ Python

.. mermaid::
    :zoom:

    flowchart TB 

        GG(["compute_grn()"])
        SS(["gen_syn_from_gf_*()"])
        EE(["compute_strain()"])
        RR(["compute_rotation()"])
        TT(["compute_stress()"])

        G["Compute Green's Functions
        (and its Spatial Derivatives)"]
        S["Compute displacements with focal mechanism
        (and its Spatial Derivatives)"]
        E["Compute Strain Tensor"]
        R["Compute Rotation Tensor"]
        T["Compute Stress Tensor"]

        GG --> G
        G --> SS --> S
        S --> EE --> E
        S --> RR --> R
        S --> TT --> T
        
        classDef cmdcls fill:#f9f2d9,stroke:#e8d174,stroke-width:2px,color:#333;
        class GG,SS,EE,RR,TT cmdcls

+ C

.. mermaid::
    :zoom:

    flowchart TB 

        GG(["grt"])
        SS(["grt.syn"])
        EE(["grt.strain"])
        RR(["grt.rotation"])
        TT(["grt.stress"])

        G["Compute Green's Functions
        (and its Spatial Derivatives)"]
        S["Compute displacements with focal mechanism
        (and its Spatial Derivatives)"]
        E["Compute Strain Tensor"]
        R["Compute Rotation Tensor"]
        T["Compute Stress Tensor"]

        GG --> G
        G --> SS --> S
        S --> EE --> E
        S --> RR --> R
        S --> TT --> T
        
        classDef cmdcls fill:#f9f2d9,stroke:#e8d174,stroke-width:2px,color:#333;
        class GG,SS,EE,RR,TT cmdcls

计算 **静态解** 的过程相同，只是函数名称和程序名称有变化。


以下是一些简易教程，可快速上手。Github主页的 :rst:dir:`example/` 文件夹中有更多示例，可在 `Github Releases <https://github.com/Dengda98/PyGRT/releases>`_ 中下载。


.. toctree::
   :maxdepth: 1

   prepare
   dynamic/gfunc
   dynamic/syn
   static/static_gfunc
   static/static_syn
   dynamic/strain_stress
   static/static_strain_stress



