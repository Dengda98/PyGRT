:author: 朱邓达
:date: 2025-12-15

波数积分的收敛性
========================

在特殊情况下，波数积分将变得难以收敛。
以下通过输出不同波数时的核函数来观察波数积分的收敛情况，
并介绍可用的积分收敛方法。

.. toctree::
   :maxdepth: 1

   integ_converg
   ptam
   dcm


可设置相应参数来指定收敛方法：

.. tabs:: 

   .. group-tab:: C 

      详见 :doc:`/Module/greenfn` 和 :doc:`/Module/static_greenfn` 模块的 **-C** 选项。

   .. group-tab:: Python

      :func:`greenfn() <pygrt.pymod.PyModel1D.greenfn>`
      函数和 :func:`static_greenfn() <pygrt.pymod.PyModel1D.static_greenfn>`
      函数支持设置 ``converg_method:str`` 来指定收敛方法（"AUTO", "NONE", "DCM", "PTAM"）。
