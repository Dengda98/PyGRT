:author: 朱邓达
:date: 2025-04-20
:updated_date: 2025-12-05

核函数 :math:`f-v` 谱图（旧版） 
==================================

在 :doc:`/Advanced/integ_converg/integ_converg` 部分介绍了程序可以导出积分过程中不同频率的核函数值，这使得我们可以绘制核函数的 :math:`f-v` 谱图 [#]_ ，以观察其中的频散特征 。

以薄层模型为例，

.. literalinclude:: run/mod1
    :language: text 

.. [#] 感谢席超强博士 `@xichaoqiang <https://github.com/xichaoqiang>`_ 的建议。


生成核函数文件
----------------------

.. note:: 

    与之前的示例相比，这里添加了控制波数积分的参数，以显式控制核函数的采样间隔和范围。但这并非是计算格林函数的最佳参数，故此时格林函数计算结果很可能是不收敛的。
    这里的目的仅是提取特定范围内不同波数下的核函数。


.. tabs:: 

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash 
            :start-after: BEGIN GRN
            :end-before: END GRN


    .. group-tab:: Python
    
        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN


读取，插值，获得 :math:`f-v` 谱
----------------------------------
Python端提供了 :py:func:`pygrt.utils.read_kernels_freqs` 函数来完成所有频率的核函数读取以及从波数 ``k`` 插值到特定速度点 ``v`` 的工作 ( :math:`v=\omega/k` )，最后获得保存有核函数的 :math:`f-v` 谱的字典。

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN read
    :end-before: END read

绘制核函数 :math:`f-v` 谱图
----------------------------------
以下将使用Python进行图件绘制。

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN plot
    :end-before: END plot

----------------------------------

+ **虚部**

.. image:: run/imag.png
    :align: center

----------------------------------

+ **实部**

.. image:: run/real.png
    :align: center


----------------------------------


从图上可看到， **核函数** :math:`f-v` **谱图的虚部峰值正是频散曲线的位置** ，尽管其中不同震源的核函数幅值不同，但 :math:`q_m，w_m` 呈现的频散特征一致，这对应的是 **Rayleigh波频散** ，而 :math:`v_m` 对应的是 **Love波频散** 。

不同的震源场点深度、不同的虚频率都会导致幅值变化，积分间隔则会影响插值效果。

.. note::

    注意其中仍然使用了虚频率，尽管这对于面波频散的分析引入了少许误差，但这使得 :math:`f-v` 谱图中的高亮区域变宽，方便观察。
