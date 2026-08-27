:author: 朱邓达
:date: 2025-09-20

计算静态解
=====================

计算 **静态解** 的主要计算流程如下：

+ Python

.. mermaid::

    flowchart TB 

        GG(["pygrt.PyModel1D.static_greenfn()"])
        SS(["pygrt.PyModel1D.static_syn()"])
        EE(["pygrt.utils.static_strain()"])
        RR(["pygrt.utils.static_rotation()"])
        TT(["pygrt.utils.static_stress()"])
        SP(["pygrt.utils.static_sproj()"])
        CO(["pygrt.utils.static_coulomb()"])

        G["Compute Green's Functions
        (and its Spatial Derivatives)"]
        S["Compute displacements with focal mechanism
        (and its Spatial Derivatives)"]
        E["Compute Strain Tensor"]
        R["Compute Rotation Tensor"]
        T["Compute Stress Tensor"]
        P["Project Stress Tensor"]
        C["Compute Coulomb Stress"]

        GG --> G
        G --> SS --> S
        S --> EE --> E
        S --> RR --> R
        S --> TT --> T
        T --> SP --> P
        P --> CO --> C
        
        classDef cmdcls fill:#FBE8CE,stroke:#BFA28C,stroke-width:2px,color:#333;
        class GG,SS,EE,RR,TT,SP,CO cmdcls

+ C (module name)

.. mermaid::

    flowchart TB 

        GG(["static_greenfn"])
        SS(["static_syn"])
        EE(["static_strain"])
        RR(["static_rotation"])
        TT(["static_stress"])
        SP(["static_sproj"])
        CO(["static_coulomb"])

        G["Compute Green's Functions
        (and its Spatial Derivatives)"]
        S["Compute displacements with focal mechanism
        (and its Spatial Derivatives)"]
        E["Compute Strain Tensor"]
        R["Compute Rotation Tensor"]
        T["Compute Stress Tensor"]
        P["Project Stress Tensor"]
        C["Compute Coulomb Stress"]

        GG --> G
        G --> SS --> S
        S --> EE --> E
        S --> RR --> R
        S --> TT --> T
        T --> SP --> P
        P --> CO --> C
        
        classDef cmdcls fill:#FBE8CE,stroke:#BFA28C,stroke-width:2px,color:#333;
        class GG,SS,EE,RR,TT,SP,CO cmdcls

.. toctree::
   :hidden:
   :maxdepth: 1

   static_gfunc
   static_syn
   static_strain_stress
   stlib
   src_fault
   rcv_fault
