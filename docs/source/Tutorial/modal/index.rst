:author: 朱邓达
:date: 2026-05-19

计算动态解（面波解）
=====================

计算 **动态面波解** 的主要计算流程如下（目前仅支持 C 模块）。其中值得注意的是，
在经过 **modsum** 模块计算得到不同震源的面波格林函数后，
后续的合成位移/应力等流程完全和 :doc:`/Tutorial/dynamic/index` 是一样的，共享计算模块。

关于面波计算的基本理论详见 |yao2026p| ，这里不再做过多介绍。

.. mermaid::

    flowchart TB 

        VV(["eigenv"])
        FN(["eigenfn"])
        MS(["modsum"])
        SS(["syn"])
        EE(["strain"])
        RR(["rotation"])
        TT(["stress"])

        V["Compute Dispersion Curves"]
        F["Compute Eigenfunctions, Energy Integrals,
           Group Velocity and Dispersion Sensitivity"]
        M["Compute Surface-wave Green's Functions (and its Spatial Derivatives)"]

        VV --> V
        V --> FN --> F
        V --> MS --> M
        M --> SS
        SS --> EE
        SS --> RR
        SS --> TT
        
        classDef cmdcls1 fill:#A1E3F9,stroke:#006BFF,stroke-width:2px,color:#333;
        classDef cmdcls2 fill:#f9f2d9,stroke:#e8d174,stroke-width:2px,color:#333;
        class VV,FN,MS cmdcls1 
        class SS,EE,RR,TT cmdcls2


.. toctree::
   :hidden:
   :maxdepth: 1

   eigenvalue
   eigenfunction
   surface_wave