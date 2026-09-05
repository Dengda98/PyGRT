:author: 朱邓达
:date: 2025-11-14

Lamb 问题
================

北京大学张海明教授已在其编著的 **《地震学中的 Lamb 问题》上下两册** 中对三类 Lamb 问题进行了详细论述，
并分别在上下两册中给出了频域解和时域解，过程之详细令人叹为观止。这里不再对 Lamb 问题及其解法做过多介绍，
详见张海明老师的书 (|zhang2021|; |zhang2024|)（强烈推荐！！）

得益于书中详细的推导过程，这里我对照下册书 **编程实现了第一、二、三类 Lamb 问题在时域的广义闭合解** 。
以下简单示范用法以及绘图，其中复现了书中的一些结果，并和频域解进行了对比。

.. toctree::
   :maxdepth: 1

   lamb1
   lamb2
   lamb3

.. note::

    **程序输入的时间为无量纲时间，输出的位移为和阶跃函数卷积后的无量纲位移。**

    + **无量纲时间** :math:`\bar{t} = \dfrac{t}{T_S} = \dfrac{t}{r/\beta} = \dfrac{\beta t}{r}`，
      其中 :math:`T_S = \dfrac{r}{\beta}` 是 S 波传播时间尺度，
      :math:`r` 为源点和接收点的两点直线距离，不是水平震中距
    + **无量纲位移** 均为 :math:`\bar{\mathbf{G}}^H = \pi^2 \mu r \mathbf{G}^H`，
      上标 :math:`H` 表示已和阶跃函数卷积，下同
    + **无量纲空间导数** 均为 :math:`\bar{\mathbf{G}}^H_{,k}=\pi^2\mu r^2\mathbf{G}^H_{,k}`

.. note::

    **如果使用了相关功能，您还需引用：**

    + Feng, X., Zhang, H., 2018. Exact closed-form solutions for lamb’s problem. Geophys. J. Int. 214, 444–459. https://doi.org/10.1093/gji/ggy131
    + Feng, X., Zhang, H., 2021. Exact closed-form solutions for lamb’s problem—III: the case for buried source and receiver. Geophys. J. Int. 224, 517–532. https://doi.org/10.1093/gji/ggaa485
    + 张海明，冯禧，2024. 地震学中的Lamb问题（下）[M]. 北京：科学出版社.
