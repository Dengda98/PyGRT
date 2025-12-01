:author: 朱邓达
:date: 2025-12-01

简介
=============

**PyGRT** 是什么？
--------------------------

**PyGRT** 是一个用于计算半无限水平层状介质中理论地震图的 **C/Python** 程序包。
**GRT** 指代的是程序计算的主要理论基础为广义反射透射系数矩阵法
(**G**\eneralized **R**\eflection-**T**\ransmission coefficient matrix)，
但不局限于该方法。

目前， **PyGRT** 可计算以下物理量（动态解和静态解）：

+ **位移全波解及其空间导数**
+ **应力、应变、旋转张量**
+ 计算 **面波** 的模块即将发布。


**PyGRT** 运行平台
---------------------------------------
**PyGRT** 发布了二进制安装包，支持在 **Windows, Linux, MacOS** 中直接下载使用，无需编译。


**PyGRT** 程序特点
---------------------------

1. **C语言的高效 + Python语言的便捷**

  + **底层复杂运算完全由 C 语言重新实现**，并基于 `OpenMP <https://www.openmp.org/>`_ 进行并行优化，极大提升计算效率。

  + C 代码被编译链接成动态库 ``libgrt.so`` ， **PyGRT** 再基于 Python 的 
    `ctypes <https://docs.python.org/3/library/ctypes.html>`_ 标准库实现对 C 库函数的调用。
    再基于第三方库 `NumPy <https://numpy.org/>`_ 、 `SciPy <https://scipy.org/>`_ 
    和 `ObsPy <https://docs.obspy.org/>`_ ，用户可很方便地完成对 C 程序结果的数据整合、
    Fourier 变换、卷积、滤波、保存到 sac 文件等操作（例如 FFT 点数不再强制要求 2 次幂）。
    借用 Python 语言的特点以及丰富成熟的第三方库，用户可灵活地实现后续的各种数据处理。

2. **C 程序 grt**

  尽管程序包名以 **Py** 打头，但仍然使用 C 语言编译了可执行文件 **grt** ，
  这使得用户既可以编写 Python 脚本，也可以在终端运行 **grt** 的形式来调用程序进行计算。

3. **模块化**
  
  + **PyGRT** 遵循模块化设计思想，不同计算模块和处理函数划分在不同模块中，
    使用不同的源文件和头文件管理。方便后续代码维护和功能扩展。
  + 受画图程序 `GMT <https://www.generic-mapping-tools.org/>`_ 的启发， 
    **grt** 程序也以模块化的形式调用不同计算模块。

4. **明确的中文注释**

  我在代码中对计算过程给出了明确的中文注释以及相关的公式索引， 
  **所有公式索引（除非特别指明）均来自**  :ref:`该初稿 <yao_init_manuscripts>` 。
  代码基本按照公式推导的计算逻辑，在不过多损失计算效率的前提下进行适当优化，保证了代码的可读性。
  若你对具体实现过程感兴趣，尤其是对方法的具体实现流程困惑时，可以在运行程序的同时结合阅读相关 **C代码** ，
  希望能给学习过程中的你提供一些参考，解答一些疑惑。仅是我浅薄的理解，仅作参考，欢迎指正。

5. **特殊源点场点分布下的计算优化** （全波解）

  - **当震中距** :math:`r` **很大时** 。由于结果中Bessel函数的形式为 :math:`J_m(kr)`，
    当震中距变大时，积分要求波数 :math:`k` 的积分间隔就越小，导致计算变慢。
    目前 **PyGRT** 实现了 **基于线性插值的Filon积分** :ref:`(纪晨, 姚振兴, 1995) <jichen_1995>` :ref:`(初稿) <yao_init_manuscripts>` 
    和 **自适应Filon积分** :ref:`(Chen and Zhang, 2001) <chen_2001>` :ref:`(张海明, 2021) <zhang_book_2021>` 以缓解此问题，可显著提高计算速度。

  - **当震源深度和台站深度很接近时（如第一类Lamb问题，台站和震源均位于地表）**。
    此时核函数随着波数 :math:`k` 的增加收敛非常慢，会在收敛值上下波动，导致很难到达到指定的收敛条件，
    需要耗费更多时间以达到更高的积分上限。
    目前 **PyGRT** 实现了以下两种方法来缓解此问题：

    + **峰谷平均法** :ref:`(Zhang et al., 2003) <zhang_2003>` :ref:`(张海明, 2021) <zhang_book_2021>` 。
      当波数 :math:`k` 达到一定上限时，可以统计波动的波峰波谷值，再递归取缩减序列 :math:`M_i \leftarrow 0.5\times(M_i + M_{i+1})` 得到收敛值。
    + **直接收敛法** :ref:`(Zhu et al., in press)` 。即将发布。

6. **开源透明，持续维护**

  **PyGRT** 开源在 `Github  <https://github.com/Dengda98/PyGRT>`_ ，将持续公开维护并扩展相关的计算功能。
  欢迎用户提交 bug 报告、申请新功能等，也欢迎提交贡献。


引用
----------------
如果你的研究中使用了 **PyGRT** 程序，请引用以下文章：

+ **PyGRT** 程序包。 由于审稿过程中程序在不断扩展完善，目前程序的功能已超过该文章的叙述范围，具体功能详见本文档。

  Zhu, D., J. Wang, J. Hao, S. Yao, Y. Xu, T. Xu, and Z. Yao (in press). 
  PyGRT: An Efficient and Integrated Python Package for Computing Synthetic 
  Seismograms in a Layered Half-Space Model, Seismological Research Letters.

+ **直接收敛法**
  
  Zhu, D., T. Xu, J. Hao, and Z. Yao (in press). A Direct Convergence Method 
  for Computing Synthetic Seismograms for a Layered Half-space with Sources and 
  Receivers at Close Depths, Bulletin of the Seismological Society of America.


主要参考
------------

相关理论方法主要参考了以下文章和书籍：

.. _bouchon_1981: 

.. [#] Michel Bouchon. 1981. A simple method to calculate Green's functions for elastic layered media. Bulletin of the Seismological Society of America. 71(4). 959–971. doi: `10.1785/BSSA0710040959 <https://doi.org/10.1785/BSSA0710040959>`_

.. _jichen_1995: 

.. [#] 纪晨, 姚振兴, 1995. 区域地震范围的宽频带理论地震图算法研究[J]. 地球物理学报, 38(4): 460-468.

.. _xie&yao_1989:

.. [#] 谢小碧, 姚振兴, 1989. 计算分层介质中位错点源静态位移场的广义反射、透射系数矩阵和离散波数方法[J]. 地球物理学报, 32(3): 270-280.

.. _yao_init_manuscripts:

.. [#] 姚振兴, 谢小碧. 2022/03. 理论地震图及其应用（初稿）.  （预计2026年出版）

.. _yao&harkrider_1983:

.. [#] Yao Z. X. and D. G. Harkrider. 1983. A generalized refelection-transmission coefficient matrix and discrete wavenumber method for synthetic seismograms. BSSA. 73(6). 1685-1699. doi: `10.1785/BSSA07306A1685 <https://doi.org/10.1785/BSSA07306A1685>`_ 

.. _zhang_book_2021: 

.. [#] 张海明 著. 2021. 地震学中的Lamb问题（上）. 科学出版社.  

.. _zhang_book_2024: 

.. [#] 张海明, 冯禧 著. 2024. 地震学中的Lamb问题（下）. 科学出版社.  

.. _zhang_2003: 

.. [#] Zhang, H. M., Chen, X. F., and Chang, S. 2003. An efficient numerical method for computing synthetic seismograms for a layered half-space with sources and receivers at close or same depths. Seismic motion, lithospheric structures, earthquake and volcanic sources: The Keiiti Aki volume, 467-486. doi: `10.1007/978-3-0348-8010-7_3 <https://doi.org/10.1007/978-3-0348-8010-7_3>`_ 

.. _chen_2001:

.. [#] Chen, X., Zhang, H., 2001. An Efficient Method for Computing Green’s Functions for a Layered Half-Space at Large Epicentral Distances. Bulletin of the Seismological Society of America 91, 858–869. doi: `10.1785/0120000113  <https://doi.org/10.1785/0120000113>`_ 

.. _kennett&kerry_1979:

.. [#] Kennett, B. L. N., and N. J. Kerry, 1979. Seismic waves in a stratified half space, Geophysical Journal International, 57, no. 3, 557–583, doi: `10.1111/j.1365-246X.1979.tb06779.x <https://doi.org/10.1111/j.1365-246X.1979.tb06779.x>`_
