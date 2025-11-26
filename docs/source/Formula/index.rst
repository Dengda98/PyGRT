公式补充
================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------

（以下面向选修中国科学院 《理论地震学》 课程的学生）

在 2022 年春季，研究生一年级的我在雁栖湖选修了中国科学院 《理论地震学》 课程，教材包括 :ref:`《理论地震图及其应用（初稿）》 <yao_init_manuscripts>` ，后因兴趣，研读了相关论文，推导公式，编写了 **PyGRT** 程序包。以下是摸索过程中推导的一些公式（持续补充），包括扩展程序中用到的，作为对教材公式的补充，欢迎批评指正。

.. note::

    公式推导过程中的Z轴取向下为正，详见 :ref:`Warning <warning_C_python_Z_direction>` 。

**反射透射系数矩阵**

.. hlist::
    :columns: 1

    - :doc:`RT`
    - :doc:`RT_static`
    - :doc:`RT_liquid`
    - :doc:`uiz`
    - :doc:`DS_zero`

**静力学震源参数**

.. hlist::
    :columns: 1

    - :doc:`static_uniform`
    - :doc:`static_source`



.. toctree::
    :maxdepth: 1
    :hidden:

    RT
    RT_static
    RT_liquid
    uiz
    DS_zero
    static_uniform
    static_source
    