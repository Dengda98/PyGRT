核函数频谱图 
=============

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn


-----------------------------------------------------------

在 :ref:`integ_converg_rst` 部分介绍了程序可以导出积分过程中不同频率的核函数值，这使得我们可以绘制核函数的频谱图 [#]_ ，以观察其中的频散特征 。

以薄层模型为例，

.. literalinclude:: run/mod1
    :language: text 

.. [#] 感谢席超强博士 `@xichaoqiang <https://github.com/xichaoqiang>`_ 的建议。


生成核函数文件
----------------------

.. note:: 

    与之前的示例相比，这里添加了控制波数积分的参数，以显式控制核函数的采样间隔和范围。但这并非是计算格林函数的最佳参数，故此时格林函数计算结果很可能是不收敛的。


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


读取，采样，获得核函数频谱
----------------------------------
Python端提供了 :py:func:`pygrt.utils.read_kernels_freqs` 函数来完成所有频率的核函数读取以及从波数 ``k`` 插值到特定速度点 ``v`` 的工作 ( :math:`k=\omega/v` )，最后获得保存有核函数频谱的字典。

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN read
    :end-before: END read

绘制核函数频谱图
------------------------
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


从图上可看到， **核函数频谱的虚部峰值正是频散曲线的位置** ，其中 :math:`q_m，w_m` 呈现的频散特征一致，只是不同震源的核函数幅值不同，这对应的是 **Rayleigh波频散** ，而 :math:`v_m` 对应的是 **Love波频散**。