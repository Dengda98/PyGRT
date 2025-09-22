计算动态全波解
=====================

计算 **动态全波解** 的主要计算流程如下：

+ Python

.. mermaid::

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

.. toctree::
   :hidden:
   :maxdepth: 1

   gfunc
   syn
   strain_stress