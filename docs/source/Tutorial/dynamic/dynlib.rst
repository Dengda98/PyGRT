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
动态格林函数按震源深度、台站深度和震中距存储，合成时只能选择库中已经存在的精确采样值，不进行插值。

**目录下要求只能有一个模型的结果，如果包含多个模型会对后续计算造成影响。**

快速上手
---------

下面的例子计算 2、4 km 两个震源深度，0、2 km 两个台站深度，以及 5、8、10 km 三个震中距，
因此会得到 12 个格林函数子目录。

.. tabs::

    .. group-tab:: C

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

        Python 接口使用列表表达相同的深度组合，并在内部调用同一个 C 模块：

        .. literalinclude:: run_library/run.py
            :language: python
            :start-after: BEGIN GRN
            :end-before: END GRN

        这里 ``grn=`` 指定的是库根目录，而不是某一个深度和距离对应的子目录。

深度和距离的写法
------------------

命令行支持三种列表形式，**-Ds**、**-Dr** 和 **-R** 的规则分别如下：

* ``-Ds2/6/2``：生成 2、4、6 km
* ``-Ds2,4,6``：直接给出逗号分隔的列表
* ``-Dsdepsrc.txt``：从文件中逐行读取

**-Dr** 和 **-R** 的写法相同。Python 中直接传入序列即可，例如
``depsrc=[2, 4, 6]``、``deprcv=[0, 2]`` 和 ``dists=[5, 8, 10]``。
程序会遍历所有震源深度和台站深度组合，再对每个组合计算所有震中距。

从库中选择格林函数
--------------------

:doc:`/Module/syn` 模块进行合成时可以把 **-G** 指向库根目录，也可以直接指向某一个子目录。
根目录模式适合脚本根据事件信息选择格林函数，此时当库中某个维度有多个值时，必须用 **-Ds**、**-Dr** 或 **-R**
明确指定；如果该维度只有一个值，可以省略对应选项。所有选择都要求与 SAC 头中的震源深度、台站深度和震中距精确相等。

.. tabs::

    .. group-tab:: C

        下面的命令从根目录中选择 ``milrow_4_2_8``，然后合成爆炸源结果：

        .. literalinclude:: run_library/run.sh
            :language: bash
            :start-after: BEGIN SYN ROOT
            :end-before: END SYN ROOT

        也可以直接指定子目录。子目录已经包含三个选择值，此时不能再设置 **-Ds/-Dr/-R**：

        .. literalinclude:: run_library/run.sh
            :language: bash
            :start-after: BEGIN SYN SUBDIR
            :end-before: END SYN SUBDIR

    .. group-tab:: Python

        Python 的 ``compute_syn()`` 使用 ``depsrc``、``deprcv`` 和 ``dist`` 传入选择值，
        其含义与 C 模块的 **-Ds/-Dr/-R** 完全一致：

        .. literalinclude:: run_library/run.py
            :language: python
            :start-after: BEGIN SYN
            :end-before: END SYN

        由于动态合成不做深度或距离插值，如果目标值不在库中，应在建库时补充该采样值，或重新计算一个单独的格林函数。
        更完整的震源类型、时间函数和分量旋转选项见 :doc:`/Module/syn`，对应的 Python 接口见
        :meth:`compute_syn() <pygrt.pymod.PyModel1D.compute_syn>`。
