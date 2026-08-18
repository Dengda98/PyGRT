:author: 朱邓达
:date: 2026-08-17

Okada 解
================

Okada 解是均匀弹性半空间中位错源引起静态变形的闭合解析解。
Okada (1985) 给出了地表位移的解析表达式，Okada (1992) 将结果推广到半空间内部，
统一处理走向滑动、倾向滑动和张裂等断层类型。NIED 提供的
`DC3D0/DC3D 程序说明 <https://www.bosai.go.jp/e/dc3d.html>`__
对点源、有限矩形断层、位移以及位移空间导数进行了整理。

**目前， Okada 解已经作为辅助模块引入到 PyGRT 中。**

.. note::

    Okada 解的适用介质是均匀各向同性弹性半空间。
    在这个特殊模型中，它作为解析解可以为 PyGRT 的数值结果提供 benchmark。


在 PyGRT 计算流程中的位置
--------------------------

PyGRT 的常规静态计算流程是先建立静态格林函数，再根据震源机制合成位移，最后由位移空间导数计算应变、旋转和应力。
对于均匀半空间，:doc:`/Module/okada` 模块提供的是解析解，
它直接从介质参数和震源几何计算位移。这相当于替代了常规流程中的
:doc:`/Module/static_greenfn` 和 :doc:`/Module/static_syn` 两个阶段，但输出仍然保留静态
合成结果的 NetCDF 接口。

因此如果还需要计算应变、旋转或应力等物理量，只需将 :doc:`/Module/okada` 生成的 NetCDF 文件交给相应的静态后处理模块即可。
换言之，:doc:`/Module/okada` 与 :doc:`/Module/static_greenfn` + :doc:`/Module/static_syn` 在流程上是平级的，
其他的后处理流程完全对 :doc:`/Module/okada` 的输出兼容。


PyGRT 如何引入 Okada 解析解
----------------------------

PyGRT 在 C 代码层面实现了 Okada 点源和矩形有限断层的计算函数，
并由 :doc:`/Module/okada` 模块负责参数解析、坐标转换和 NetCDF 输出。
计算过程可以概括为：

#. 将 PyGRT 的北向、东向、深度坐标转换为 Okada 局部坐标。Okada 的局部 ``X`` 沿断层
   走向，``Y`` 为上倾方向的水平投影，``Z`` 向上
#. 将爆炸源、双力偶源和张裂源转换为 Okada 的 potency；将 Coulomb 断层记录转换为
   矩形断层的长度、宽度、倾角、滑动和顶深
#. 调用点源解或有限断层解，得到局部坐标下的位移及位移偏导
#. 将结果旋转回 PyGRT 的 ZRT 或 ZNE 分量，并按静态模块约定换算为 cm 和无量纲导数
#. 将结果写入与 :doc:`/Module/static_syn` 相同风格的 NetCDF 文件


命令行与 Python 接口
----------------------

命令行模块的参数说明见 :doc:`/Module/okada`。Python 中可以使用
:func:`compute_okada() <pygrt.utils.compute_okada>`，其输入参数与
:func:`compute_static_syn() <pygrt.pymod.PyModel1D.compute_static_syn>` 基本一致，
不过需指定半空间模型的 P 波速度、S 波速度和密度。
其中点源类型由 ``strike``、``dip`` 和 ``rake`` 自动确定：三个参数均不设置时为爆炸源，
设置 ``strike`` 和 ``dip`` 时为张裂源，同时设置 ``rake`` 时为剪切源。


输入与输出的兼容性
--------------------

:doc:`/Module/okada` 在输入组织、震源和接收点表达方式上与 :doc:`/Module/static_syn`
保持对齐，仅将静态格林函数输入替换为均匀半空间的 *vp*、*vs* 和 *rho* 参数，
并直接计算点源或有限断层的静态位移。

输出同样采用与 :doc:`/Module/static_syn` 对齐的 |NetCDF| 格式，包括位移、可选的位移空间导数以及坐标和介质属性。
因此，Okada 结果可以直接应用到后续静态计算流程。


示例
-------

下面分别计算一个点源和一条有限矩形断层的位移场。示例使用半空间
:math:`V_P=6.0` km/s、:math:`V_S=3.464` km/s、:math:`\rho=2.7` g/cm\ :sup:`3`，
并将结果输出为 ZNE 分量。

点源命令如下：

.. literalinclude:: run/run.sh
    :language: bash
    :start-after: BEGIN POINT
    :end-before: END POINT

有限断层命令如下：

.. literalinclude:: run/run.sh
    :language: bash
    :start-after: BEGIN FINITE
    :end-before: END FINITE

点源和有限断层的位移结果分别为：

.. figure:: run/okada_point.svg
    :align: center

    点源 Okada 解的 Z、N、E 位移

.. figure:: run/okada_finite.svg
    :align: center

    有限矩形断层 Okada 解的 Z、N、E 位移

完整脚本见 :download:`run.sh <run/run.sh>` 和 :download:`plot.py <run/plot.py>`。


参考文献
----------

* Okada, Y. (1985). Surface deformation due to shear and tensile faults in a half-space.
  *Bulletin of the Seismological Society of America*, 75, 1135–1154
* Okada, Y. (1992). Internal deformation due to shear and tensile faults in a half-space.
  *Bulletin of the Seismological Society of America*, 82, 1018–1040
