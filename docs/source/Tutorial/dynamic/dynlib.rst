:author: 朱邓达
:date: 2026-08-17

创建动态全波格林函数库
========================================

动态格林函数可以重复用于许多不同的震源机制和方位角。实际计算时，模型、震源深度、台站深度和震中距确定后，
格林函数就确定了，因此在模型固定后，可以把常用的深度和距离一次性计算成一个库。这里的“创建库”不需要再手动编写循环，
而是把深度列表和震中距列表传给 :doc:`/Module/greenfn` 模块，程序在内部遍历所有组合。

动态格林函数使用 SAC 文件保存，目录结构为
``{outdir}/{model}_{depsrc}_{deprcv}_{dist}/{stype}.sac``。
例如 ``GRN/milrow_4_2_8/`` 表示震源深度 4 km、台站深度 2 km、震中距 8 km。
动态格林函数按震源深度、台站深度和震中距存储，合成时默认对库内深度和距离插值，也可选择最近邻查询。

**目录下要求只能有一个模型的结果，必须包含震源深度、接收深度和震中距三个采样轴的全部组合，
各轴可以非等距。**
如需合成空间导数，建库时须使用 **-e** 或 *calc_upar=True*。

快速上手
---------

下面的例子计算 2、4 km 两个震源深度，0、2 km 两个台站深度，以及 5、8、10 km 三个震中距，
因此会得到 12 个格林函数子目录。

.. tabs::

    .. group-tab:: CLI

        .. literalinclude:: run_library/run.sh
            :language: bash
            :start-after: BEGIN GRN
            :end-before: END GRN

        运行后可以看到类似下面的目录：

        .. code-block:: text

            GRN/
            ├── milrow_2_0_5/
            ├── milrow_2_0_8/
            ├── milrow_2_0_10/
            ├── milrow_4_0_5/
            └── ...

        **-Ds** 和 **-Dr** 必须成对使用。深度列表和 **-R** 的距离列表都必须严格递增，
        也可以使用等距范围或每行一个数值的文件，具体语法见 :doc:`/Module/greenfn`。

    .. group-tab:: Python

        Python 接口使用列表表达相同的深度组合，并在内部调用同一个 :command:`grt` 模块：

        .. literalinclude:: run_library/run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

        这里 ``grn=`` 指定的是库根目录，而不是某一个深度和距离对应的子目录。

深度和距离的写法
------------------

**-Ds**、**-Dr** 和 **-R** 命令行支持三种列表形式，规则如下：

* ``-Ds2/6/2``：生成 2、4、6 km
* ``-Ds2,4,6``：直接给出逗号分隔的列表
* ``-Dsdepsrc.txt``：从文件中逐行读取

**-Dr** 和 **-R** 的写法相同。Python 中直接传入序列即可，例如
``depsrc=[2, 4, 6]``、``deprcv=[0, 2]`` 和 ``dists=[5, 8, 10]``。
程序会遍历所有震源深度和台站深度组合，再对每个组合计算所有震中距。

从格林函数库合成
--------------------

使用 :doc:`/Module/syn` 合成时，将 **-G** 指向库根目录，指定目标源深、台深和震中距即可。
目标位置不必与库中的采样节点重合，只需位于库的采样范围内。
可选择以下两种方式：

* **线性插值（默认）**：根据相邻节点的合成结果插值，CLI 使用 **-i1**，Python 使用 *interpolate=True*
* **最近邻**：使用最近的库节点，CLI 使用 **-i0**，Python 使用 *interpolate=False*

若 **-G** 直接指定单个节点子目录，则直接使用该节点，无需选择插值方式。

例如，上面建立的库覆盖源深 2～4 km、台深 0～2 km、震中距 5～10 km。
下面在源深 3 km、台深 1 km、震中距 7 km 处合成爆炸源结果，
该位置并非库中的节点，两种查询方式均可使用：

.. tabs::

    .. group-tab:: CLI

        **-Ds/-Dr/-R** 分别指定目标源深、台深和震中距，两次合成仅改变查询方式和输出目录：

        .. literalinclude:: run_library/run.sh
            :language: bash
            :start-after: BEGIN SYN
            :end-before: END SYN

    .. group-tab:: Python

        *depsrc/deprcv/dist* 分别指定目标源深、台深和震中距，用 *interpolate* 选择查询方式：

        .. literalinclude:: run_library/run.py
            :language: python
            :start-after: BEGIN SYN
            :end-before: END SYN

