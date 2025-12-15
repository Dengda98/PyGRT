:author: 朱邓达
:date: 2025-04-19

峰谷平均法
========================

具体原理很简洁易懂，从以下图像中你也能了解大概。详见 :ref:`(Zhang et al., 2003) <zhang_2003>` :ref:`(张海明, 2021) <zhang_book_2021>`。

在 :doc:`integ_converg` 介绍中，当震源和场点深度接近时，你会注意到输出的核函数文件增加了，
这是因为程序内部设定 **当震源和场点深度差小于1km时，自动使用峰谷平均法（PTAM）**。
具体流程为，程序中在进行完离散波数积分后，继续增加k（不同震中距使用不同dk），
使用PTAM寻找足够数量的波峰和波谷（内部设定为36个），
再对这些波峰波谷取缩减序列 :math:`M_i \leftarrow 0.5\times(M_i + M_{i+1})` ，得到估计的积分收敛值。
取缩减序列的C代码如下，

.. code-block:: C 

    for(int n=n1; n>1; --n){
        for(int i=0; i<n-1; ++i){
            arr[i] = 0.5*(arr[i] + arr[i+1]);
        }
    }


在 ``K_{iw}_{freq}`` 文件同级目录下，程序把 **PTAM过程中的核函数以及积分峰谷位置分为两个文件** 保存在 ``PTAM_{ir}_{dist}/`` 目录下（ ``{ir}`` 为震中距索引， ``{dist}`` 为震中距），其中  ``PTAM_{ir}_{dist}/K_{iw}_{freq}`` 为核函数文件（格式不变）， ``PTAM_{ir}_{dist}/PTAM_{iw}_{freq}`` 为峰谷文件，其中记录积分值的峰谷。

.. note:: 

    :doc:`/Module/ker2asc` 模块也支持将 ``PTAM_{ir}_{dist}/PTAM_{iw}_{freq}`` 文件转为文本格式，

    .. literalinclude:: run/run.sh
        :language: bash
        :start-after: BEGIN grt.ker2asc ptam
        :end-before: END grt.ker2asc ptam

    输出的文件如下，

    .. literalinclude:: run/ptam_stats_head
        :language: text

    记录了不同震源、不同积分类型的峰谷位置。


故为了绘制完整的波数积分+PTAM过程，涉及三个文件。不过Python函数已经做了优化，由于程序输出的目录结构固定，文件命名格式固定，只需传入 ``PTAM_{ir}_{dist}/PTAM_{iw}_{freq}`` 文件，Python函数会自动读取对应的三个文件。

.. literalinclude:: run/run.py
    :language: python
    :start-after: BEGIN plot ptam
    :end-before: END plot ptam

.. image:: run/SS_0_0.0_ptam_RI.png
    :align: center


红十字为选取的波峰波谷，再取缩减序列即可得到估计的积分收敛值。


静态解的积分收敛性
-------------------------
以上部分是以动态解为例，静态解从积分类型、收敛特征、文件格式、绘图完全类似，只是不再有频率索引值。

.. note:: 

    核函数文件中记录的值非最终核函数值。对于静态解，还需乘 :math:`\left(\dfrac{\Delta k}{4\pi\mu}\right)`。

假设震源深度0.1km，场点位于地表，场点仅定义一个点(2,2)做示例，这里直接给出脚本。

.. tabs::  

    .. group-tab:: C 

        .. literalinclude:: run/run.sh
            :language: bash
            :start-after: BEGIN SGRN
            :end-before: END SGRN

    .. group-tab:: Python

        .. literalinclude:: run/run.py
            :language: python
            :start-after: BEGIN SGRN
            :end-before: END SGRN

----------------------------------


+ **只使用离散波数积分**

.. image:: run/SS_0_0.1_static.png
    :align: center

----------------------------------


+ **使用峰谷平均法**
  
.. image:: run/SS_0_0.1_ptam_static.png
    :align: center


