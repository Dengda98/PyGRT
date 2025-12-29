:author: 朱邓达
:date: 2025-04-20

核函数 :math:`f-v` 谱图 
=========================

在 :doc:`/Advanced/integ_converg/integ_converg` 部分介绍了程序可以导出积分过程中不同频率的核函数值，
这使得我们可以绘制核函数的 :math:`f-v` 谱图 [#]_ ，以观察其中的频散特征 。

在本页文档的旧版本中（ :doc:`/Advanced/kernel_old/kernel` ），我们使用了格林函数计算过程中导出的核函数文件，
其中每个采样点是在波数上等距采样，需要经过插值才能转到相速度，且定义不同频率、不同相速度的采样点不方便。

为了更方便地导出核函数文件，在 0.14.0 版本起，我新增了模块 :doc:`/Module/kernel` 直接面向这一任务。

以薄层模型为例，

.. literalinclude:: run/mod1
    :language: text 

.. [#] 感谢席超强博士 `@xichaoqiang <https://github.com/xichaoqiang>`_ 的建议。


核函数文件
----------------------

使用 :command:`grt` 程序调用 :doc:`/Module/kernel` 计算核函数，

.. literalinclude:: run/run.sh
    :language: bash 
    :start-after: BEGIN KERNEL
    :end-before: END KERNEL

得到的每个频率的核函数文件名开头为 ``C`` ，表示采样位置直接记录的相速度，
且 :doc:`/Advanced/integ_converg/integ_converg` 介绍的核函数文件的读取方法依然适用，
只是对应波数的列会改为等距的相速度，列名用 “c” 表示。

.. note:: 

    :doc:`/Module/kernel` 模块保存的核函数是理论核函数值，
    这与 :doc:`/Advanced/integ_converg/integ_converg` 中介绍的在计算格林函数过程中保存的核函数有幅值差异。


.. note::

    注意其中使用了虚频率，这使得 :math:`f-v` 谱图中的高亮区域变宽，方便观察。



绘制核函数 :math:`f-v` 谱图
----------------------------------
以下将使用 Python 进行图件绘制。

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN plot
    :end-before: END plot

----------------------------------

+ **虚部**

.. figure:: run/imag.svg
    :align: center

----------------------------------

+ **实部**

.. figure:: run/real.svg
    :align: center


----------------------------------


从图上可看到， **核函数** :math:`f-v` **谱图的虚部峰值正是频散曲线的位置** ，尽管其中不同震源的核函数幅值不同，但 :math:`q_m，w_m` 呈现的频散特征一致，这对应的是 **Rayleigh波频散** ，而 :math:`v_m` 对应的是 **Love波频散** 。

用户还可调整震源场点深度观察不同部分的能量变化，还可调整相速度范围观察泄露模 （leaky-mode）。
程序包为理论研究提供了便捷工具。
