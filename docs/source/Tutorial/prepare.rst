准备工作
=============

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

**PyGRT** 程序包由C和Python两个编程语言的代码组成，目的是兼并高效性和便捷性。
底层复杂运算由C语言编写，编译链接成动态库 ``libgrt.so`` 供Python调用。Python通过 ``ctypes`` 库导入动态库以使用外部函数，以此兼并了C语言的高效和Python语言的便捷。

除了Python脚本式运行， **PyGRT** 保留传统命令行式运行C程序 :command:`grt` 。
受 `GMT <https://www.generic-mapping-tools.org/>`_ 的启发， :command:`grt` 程序对多个计算功能使用模块化管理，
可使用类似于以下格式来运行程序。每个模块可使用 ``-h`` 查看帮助。

.. code-block:: C

    grt <module-name> [<module-options>]

**C程序的运行独立于Python，不需要Python环境，从而满足了更多计算场景。** 

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

+ C (module name)

.. mermaid::
    :zoom:

    flowchart TB 

        GG(["greenfn"])
        SS(["syn"])
        EE(["strain"])
        RR(["rotation"])
        TT(["stress"])

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

计算 **静态解** 的过程相同，只是在名称上有如下区别：

+ Python 函数名称增加 "static\_" ，例如 :func:`compute_static_grn() <pygrt.pymod.PyModel1D.compute_static_grn>` 。
+ C 程序模块名增加 "static\_" 前缀，例如 :command:`static_greenfn` ；或者在 :command:`grt` 与无 "static\_" 前缀的模块名之间增加 "static" 命令，
  即支持以下两种方式（以 :command:`greenfn` 模块为例）：

  .. code-block:: bash

    grt static_greenfn [<module-options>]
    grt static greenfn [<module-options>]

-----------------

入门部分包含一些示例与说明，可快速上手。Github主页的 :rst:dir:`example/` 文件夹中有更多示例，可在 `Github Releases <https://github.com/Dengda98/PyGRT/releases>`_ 中下载。


建立模型文件
--------------------

.. image:: mod.svg
   :align: center

这里的选项卡 **C** 不代表内部是C语言代码，而是 **使用C程序** :command:`grt` 。后续的选项卡沿用此设定，不再解释。

.. tabs:: 

    .. tab:: C

        **PyGRT** 以如下自由格式定义模型中每层的物性参数，每列之间以空格隔开（最后两列的 Qp, Qs 可省略）

        .. code-block:: text

            Thickness(km)    Vp(km/s)    Vs(km/s)   Rho(g/cm^3)   [Qp]   [Qs]

        
        例如 :file:`milrow` 模型（假设文本文件名为 `milrow` ）

        .. literalinclude:: dynamic/run/milrow
            :language: text
        

    .. tab:: Python

        模型格式与C一致，在Python中可以使用 :code:`np.loadtxt()` 导入文本文件，或者手动定义数组

        .. literalinclude:: dynamic/run/run.py
            :language: python
            :start-after: START BUILD MODEL
            :end-before: END BUILD MODEL


.. note::

    最后一行表示半空间，对应厚度值不会被使用。

.. note::

    Vs 设置为 0 表示该层为液体层。