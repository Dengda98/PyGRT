:author: 朱邓达
:date: 2025-05-10

公式补充
================

（以下面向选修中国科学院 《理论地震学》 课程的学生）

在 2022 年春季，研究生一年级的我在雁栖湖选修了中国科学院 《理论地震学》 课程，使用的教材包括 |yao2026p| 的初稿，
后因兴趣，研读了相关论文，推导公式，编写了 **PyGRT** 程序包。以下是摸索过程中推导的一些公式（持续补充），包括一些扩展部分的公式，作为对教材公式的补充，欢迎批评指正。

.. note::

    公式推导过程中的Z轴取向下为正，详见 :ref:`Warning <warning_C_python_Z_direction>` 。

**反射透射系数矩阵**

.. hlist::
    :columns: 1

    - :doc:`RT/RT`
    - :doc:`RT/RT_static`
    - :doc:`RT/RT_liquid`
    - :doc:`RT/boundary`

**震源系数**

.. hlist::
    :columns: 1

    - :doc:`source/static_uniform`
    - :doc:`source/static_source`

**波数积分**

.. hlist::
    :columns: 1

    - :doc:`k-integral/DCM_supply`

**其它**

.. hlist::
    :columns: 1

    - :doc:`others/uiz`
    - :doc:`others/DS_zero`



.. toctree::
    :maxdepth: 1
    :hidden:

    RT/RT
    RT/RT_static
    RT/RT_liquid
    RT/boundary
    
    source/static_uniform
    source/static_source

    k-integral/DCM_supply`

    others/uiz
    others/DS_zero
    