:author: 朱邓达
:date: 2025-09-22


模块手册
##################

以下是 C 程序 :command:`grt` 的模块手册，包含模块各个选项的使用说明以及示例。

**动态解**

+ *全波解* （波数积分, Wavenumber Integration）

    .. mermaid::

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

    .. hlist::
        :columns: 4

        - :doc:`greenfn`
        - :doc:`syn`
        - :doc:`strain`
        - :doc:`stress`
        - :doc:`rotation`
        - :doc:`kernel`

+ *面波解* （模态叠加, Modal Summation）

    release later.


**静态解**

.. mermaid::

    flowchart TB 

        GG(["static_greenfn"])
        SS(["static_syn"])
        EE(["static_strain"])
        RR(["static_rotation"])
        TT(["static_stress"])

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

.. hlist::
    :columns: 4

    - :doc:`static_greenfn`
    - :doc:`static_syn`
    - :doc:`static_strain`
    - :doc:`static_stress`
    - :doc:`static_rotation`


**辅助模块**

.. hlist::
    :columns: 4

    - :doc:`sac2asc`
    - :doc:`ker2asc`
    - :doc:`travt`
    - :doc:`lamb1`


.. toctree::
    :maxdepth: 1
    :hidden:

    greenfn
    syn
    strain
    stress
    rotation
    kernel
    static_greenfn
    static_syn
    static_strain
    static_stress
    static_rotation
    sac2asc
    ker2asc
    travt
    lamb1

