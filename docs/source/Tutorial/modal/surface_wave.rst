:author: 朱邓达
:date: 2026-05-20

计算面波格林函数
============================================

根据留数定理，本征值处对波数积分的贡献可以解析地表达。将多个本征值得到的贡献叠加起来就得到面波格林函数，
这个方法也被称为 **模态叠加法（Modal Summation Method）** 。 |yao2026p| 给出了推导过程和具体公式，这里不再重复。

C 模块为 :doc:`/Module/modsum` （Python函数将于后续版本增加），
计算将基于 :doc:`/Module/eigenv` 模块的频散结果以及模型路径。
:doc:`/Module/modsum` 模块保存的波形格式与 :doc:`/Module/greenfn` 模块一致。

示例程序
-----------

假设在 :file:`milrow` 模型中，震源深度 2 km，接收点位于地表，震中距为 100 km，
对比位错源的全波解和不同阶面波解的波形。
其中面波波形的时域长度和间隔通过 :doc:`/Module/eigenv` 模块中频域的设置间接确定，
即 *T = 1 / df* , *dt = 1 / (2\*fmax)* 。

.. literalinclude:: run_modsum/run.sh
    :language: bash
    :start-after: BEGIN
    :end-before: END

:download:`compare_wave.py <run_modsum/compare_wave.py>` （虚线为全波解，实线为面波解）

+ **0 阶面波**

.. figure:: run_modsum/compare_0.svg

+ **1 阶面波**

.. figure:: run_modsum/compare_1.svg

+ **2 阶面波**

.. figure:: run_modsum/compare_2.svg


+ **所有阶面波**

.. figure:: run_modsum/compare_2.svg