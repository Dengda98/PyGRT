:author: 朱邓达
:date: 2026-05-20

计算频散曲线
====================

频散曲线上每个点对应的相速度被称为面波本征值。具体地，这里我们使用自适应搜根方法 (|ars2026|) 计算面波相速度频散曲线，
C 模块为 :doc:`/Module/eigenv` （Python函数将于后续版本增加）。

以下示例使用如下层状模型，文件命名为 :file:`mod2` :

.. literalinclude:: run/mod2
    :language: text


示例程序
-----------

计算 :file:`mod2` 模型中的 Rayleigh 波和 Love 波的所有阶频散曲线，频率 5 Hz 以下，间隔 0.01 Hz。

.. literalinclude:: run/run.sh
    :language: bash
    :start-after: BEGIN DISPERSION
    :end-before: END DISPERSION

计算结果分别导出为 |NetCDF| 经典格式，其中具体变量详见 :doc:`/Module/eigenv` 模块文档。

Python 提供了 :func:`read_dispersion() <pygrt.utils.read_dispersion>`
将结果读入字典（``freq``、``c`` 等按阶号分组），以及 :func:`plot_dispersion() <pygrt.utils.plot_dispersion>` 进行简易绘图。
如果仅想要快速浏览，也可使用 :doc:`/Module/disp2asc` 模块将 |NetCDF|
格式的频散转为文本格式，方便简单阅读结果以及使用 |GMT| 绘制。

以下提供这两种方式的简易绘图脚本。

:download:`plot_dispersion.py (基于 Matplotlib) <run/plot_dispersion.py>` | 
:download:`plot_dispersion.sh (基于 GMT) <run/plot_dispersion.sh>`

.. figure:: run/dispersion_py.svg

绘制久期函数
---------------------

计算本征值的过程可以转变为搜索久期函数零点的过程 (|yao2026|)。
如果感兴趣，可以使用 :doc:`/Module/eigenv` 模块导出搜根过程中对应久期函数幅值。
例如，以上述程序为例，导出频率 1 Hz 处 Rayleigh 波和 Love 波的久期函数，相速度范围在 [3.0, 5.5] km/s :

.. literalinclude:: run/run.sh
    :language: bash
    :start-after: BEGIN SECFUNC
    :end-before: END SECFUNC

Python 提供了 :func:`read_secfunc() <pygrt.utils.read_secfunc>`
将结果读入字典列表（每一条对应一层的久期函数），以及
:func:`plot_secfunc() <pygrt.utils.plot_secfunc>` 绘制其中一条（也可传入 ``{标题: 字典}`` 分栏绘制）。

:download:`plot_secfunc.py <run/plot_secfunc.py>`

.. figure:: run/secfunc.svg