:author: 朱邓达
:date: 2025-09-22


模块手册
##################

以下是 C 程序 :command:`grt` 的模块手册，包含模块各个选项的使用说明以及示例。

**动态解**

+ *全波解* （波数积分, Wavenumber Integration）

.. grid:: 3
    :gutter: 0

    .. grid-item::
        :columns: 3

    .. grid-item::
        :columns: 5

        .. container:: small-mermaid

            .. mermaid::

                flowchart LR 

                GG(["greenfn"])
                SS(["syn"])
                EE(["strain"])
                RR(["rotation"])
                TT(["stress"])

                GG --> SS
                SS --> EE
                SS --> RR
                SS --> TT
                
                classDef cmdcls fill:#f9f2d9,stroke:#e8d174,stroke-width:2px,color:#333;
                class GG,SS,EE,RR,TT cmdcls

    .. grid-item::
        :columns: 2
        :child-direction: column
        :child-align: center

        - :doc:`greenfn`
        - :doc:`syn`
        - :doc:`strain`
        - :doc:`stress`
        - :doc:`rotation`
        - :doc:`kernel`
    

+ *面波解* （模态叠加, Modal Summation）

.. grid:: 3
    :gutter: 0
    
    .. grid-item::
        :columns: 2

    .. grid-item::
        :columns: 7

        .. container:: small-mermaid

            .. mermaid::

                flowchart LR 

                VV(["eigenv"])
                FN(["eigenfn"])
                MS(["modsum"])
                SS(["syn"])
                EE(["strain"])
                RR(["rotation"])
                TT(["stress"])

                VV --> FN
                VV --> MS
                MS --> SS
                SS --> EE
                SS --> RR
                SS --> TT
                
                classDef cmdcls1 fill:#A1E3F9,stroke:#006BFF,stroke-width:2px,color:#333;
                classDef cmdcls2 fill:#f9f2d9,stroke:#e8d174,stroke-width:2px,color:#333;
                class VV,FN,MS cmdcls1 
                class SS,EE,RR,TT cmdcls2

    .. grid-item::
        :columns: 2
        :child-direction: column
        :child-align: center

        - :doc:`eigenv`
        - :doc:`eigenfn`
        - :doc:`modsum`


**静态解**

.. grid:: 3
    :gutter: 0

    .. grid-item::
        :columns: 2

    .. grid-item::
        :columns: 7

        .. container:: small-mermaid

            .. mermaid::

                flowchart LR 

                GG(["static_greenfn"])
                SS(["static_syn"])
                EE(["static_strain"])
                RR(["static_rotation"])
                TT(["static_stress"])

                GG --> SS
                SS --> EE
                SS --> RR
                SS --> TT
                
                classDef cmdcls fill:#FBE8CE,stroke:#BFA28C,stroke-width:2px,color:#333;
                class GG,SS,EE,RR,TT cmdcls

    .. grid-item::
        :columns: 2
        :child-direction: column
        :child-align: center

        - :doc:`static_greenfn`
        - :doc:`static_syn`
        - :doc:`static_strain`
        - :doc:`static_stress`
        - :doc:`static_rotation`


**辅助模块**

.. hlist::
    :columns: 4

    - :doc:`disp2asc`
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
    eigenv
    eigenfn
    modsum
    static_greenfn
    static_syn
    static_strain
    static_stress
    static_rotation
    disp2asc
    sac2asc
    ker2asc
    travt
    lamb1

