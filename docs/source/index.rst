**PyGRT** 文档
=====================

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------


.. raw:: html 

  <!-- Place this tag in your head or just before your close body tag. -->
  <script async defer src="https://buttons.github.io/buttons.js"></script>

  <!-- Place this tag where you want the button to render. -->
  <a class="github-button" href="https://github.com/Dengda98/PyGRT" data-icon="octicon-star" data-size="large" aria-label="Star Dengda98/PyGRT on GitHub">Star</a>



**PyGRT** : 一个用于计算在半无限水平分层均匀介质模型中的理论地震图的 **Python**\程序包，使用广义反射透射系数矩阵法(**G**\eneralized **R**\eflection-**T**\ransmission coefficient matrix) 以及离散波数法等方法计算。   

**PyGRT** 能计算什么？
--------------------------
目前，**PyGRT** 可计算区域尺度内半无限水平分层均匀介质中的 **动态与静态的位移，以及应力、应变、旋转张量。**


运行平台
---------
**Windows, Linux, MacOS**


程序特点
---------

+ **格林函数频谱由C程序并行计算**。由于不同频率间的计算相互独立，在计算过程中基于 `OpenMP <https://www.openmp.org/>`_ 进行了 **并行优化**，极大提升了计算效率；
\

+ **将Python语言的便携、可扩展性与C语言的计算高效特点结合**。C代码被编译链接成动态库 ``libgrt.so`` ， **PyGRT** 再基于Python的 `ctypes <https://docs.python.org/3/library/ctypes.html>`_ 标准库实现对C库函数的调用。再基于第三方库 `NumPy <https://numpy.org/>`_ 、 `SciPy <https://scipy.org/>`_ 和 `ObsPy <https://docs.obspy.org/>`_ ，用户可很方便地完成对C程序结果的数据整合、Fourier变换、卷积、滤波、保存到sac文件等操作（例如FFT点数不再强制要求2次幂）。借用Python语言的特点以及丰富成熟的第三方库，可灵活地实现后续的各种数据处理。
\

+ **明确的中文注释**。我在代码中对计算过程给出了明确的中文注释以及相关的公式索引， **所有公式索引（除非特别指明）均来自**  :ref:`该初稿 <yao_init_manuscripts>`，在代码中基本按照公式推导的计算逻辑，在不过多损失计算效率的前提下进行适当优化，保证了代码的可读性，希望能给学习过程中的你提供一些参考，解答一些疑惑。仅是我浅薄的理解，仅作参考，欢迎指正。
\

+ **特殊情况下的计算优化**。目前包括两个部分：  
\

  - **当震中距** :math:`r` **很大时**。由于结果中Bessel函数的形式为 :math:`J_m(kr)`，当震中距变大时，积分要求波数 :math:`k` 的积分间隔就越小，导致计算变慢。目前 **PyGRT** 实现了 **基于线性插值的Filon积分** :ref:`(纪晨, 姚振兴, 1995) <jichen_1995>` :ref:`(初稿) <yao_init_manuscripts>` 和 **自适应Filon积分** :ref:`(Chen and Zhang, 2001) <chen_2001>` :ref:`(张海明, 2021) <zhang_book_2021>` 以缓解此问题，可显著提高计算速度。
\

  - **当震源深度和台站深度很接近时（如第一类Lamb问题，台站和震源均位于地表）**。此时核函数随着波数 :math:`k` 的增加收敛非常慢，会在收敛值上下波动，导致很难到达到指定的收敛条件，需要耗费更多时间以达到更高的积分上限。目前 **PyGRT** 实现了 **峰谷平均法** 以缓解此问题 :ref:`(Zhang et al., 2003) <zhang_2003>` :ref:`(张海明, 2021) <zhang_book_2021>`，当波数 :math:`k` 达到一定上限时，可以统计波动的波峰波谷值，再递归取缩减序列 :math:`M_i \leftarrow 0.5\times(M_i + M_{i+1})` 得到收敛值。
\

+ **积分过程可输出记录文件**。通过函数参数设置，在计算格林函数频谱时可输出不同频率的积分过程文件，记录了核函数随波数 :math:`k` 的变化，且 **PyGRT** 包含了读取该文件的Python函数，可以简易地观察到相关量的变化，方便各种测试以及学习过程中的理解。  


--------------------------------------------------------------------------


以下包括 **Python** 和 **C** 两部分的API介绍。如果只关心程序的使用，可以只查看 **Python** 部分；若对具体实现过程感兴趣，尤其是在学习 **广义反射透射系数矩阵、离散波数法、线性Filon积分法、峰谷平均法** 等内容感到迷惑难解时，可以在运行程序的同时结合阅读相关 **C代码**，我在C代码中给出了 **明确的中文注释以及相关公式索引** ，希望能帮助你在学习相关内容时解答一些疑惑。

本文档页面源码位于 Github主页的 `docs/ <https://github.com/Dengda98/PyGRT/tree/main/docs/source>`_，可查看示例代码和画图脚本。



.. toctree::
  :maxdepth: 1
  :caption: 目录
  :numbered: 3


  文档首页 <self>
  install
  Tutorial/index 
  Advanced/index 
  API/api
  copyright

--------------------------------------------------------------------------


主要参考
---------

.. _bouchon_1981: 

.. [1] Michel Bouchon. 1981. A simple method to calculate Green's functions for elastic layered media. Bulletin of the Seismological Society of America. 71(4). 959–971. doi: `10.1785/BSSA0710040959 <https://doi.org/10.1785/BSSA0710040959>`_

.. _jichen_1995: 

.. [2] 纪晨, 姚振兴, 1995. 区域地震范围的宽频带理论地震图算法研究[J]. 地球物理学报, 38(4): 460-468.

.. _xie&yao_1989:

.. [3] 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、透射系数矩阵和离散波数方法[J]. 地球物理学报, 32(3): 270-280.

.. _yao_init_manuscripts:

.. [4] 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  

.. _yao&harkrider_1983:

.. [5] Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699. doi: `10.1785/BSSA07306A1685 <https://doi.org/10.1785/BSSA07306A1685>`_ 

.. _zhang_book_2021: 

.. [6] 张海明. 2021. 地震学中的Lamb问题（上）. 科学出版社.  

.. _zhang_2003: 

.. [7] Zhang, H. M., Chen, X. F., and Chang, S. 2003. An efficient numerical method for computing synthetic seismograms for a layered half-space with sources and receivers at close or same depths. Seismic motion, lithospheric structures, earthquake and volcanic sources: The Keiiti Aki volume, 467-486. doi: `10.1007/978-3-0348-8010-7_3 <https://doi.org/10.1007/978-3-0348-8010-7_3>`_ 

.. _chen_2001:

.. [8] Chen, X., Zhang, H., 2001. An Efficient Method for Computing Green’s Functions for a Layered Half-Space at Large Epicentral Distances. Bulletin of the Seismological Society of America 91, 858–869. doi: `10.1785/0120000113  <https://doi.org/10.1785/0120000113>`_ 