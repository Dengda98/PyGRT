
.. _install_section:

安装
=============

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------


.. |gr| replace:: `Github Releases <https://github.com/Dengda98/PyGRT/releases>`__

依赖
------------

+ 如果你想使用C程序，你需要 `Seismic Analysis Code (SAC) <http://www.iris.edu/ds/nodes/dmc/forms/sac/>`_ ，因为C程序以二进制的SAC格式输出波形。
+ 如果你想以Python脚本形式使用，依赖已在 :file:`setup.py` 中写好，直接使用 :command:`pip` 安装即可。


免编译安装
------------

目前 **PyGRT** 已分发内置预编译二进制文件的安装包，用户仅需运行以下命令即可

.. code-block:: bash

    pip install pygrt-kit

或者从 |gr| 中下载最新版本的程序压缩包（符合自己的操作系统），解压，在根目录下运行以下命令即可

.. code-block:: bash

    pip install .  

如果你不想使用Python，只想使用传统的命令行形式运行C程序，也可从 |gr| 中下载最新版本的程序压缩包（符合自己的操作系统），解压，其中 :rst:dir:`pygrt/C_extension/bin` 和 :rst:dir:`pygrt/C_extension/lib` 为预编译好的可执行文件目录和动态/静态库目录，按自己习惯配置环境变量 :envvar:`PATH` 即可。



环境变量配置
-------------
如果你使用 :command:`pip` 安装后，想使用其中的C程序（ :command:`grt` , :command:`grt.syn` etc. ），需配置环境变量 :envvar:`PATH` 。运行以下命令

.. code-block:: bash

    python -m pygrt.print

输出

.. code-block:: text
    
    PyGRT installation directory: </path/to/installation>
    PyGRT executable file directory: </path/to/installation/bin>
    PyGRT library directory: </path/to/installation/lib>

将其中的“PyGRT executable file directory”路径添加到环境变量 :envvar:`PATH` 中即可。C程序的运行独立于Python，每个C程序可使用 ``-h`` 查看帮助。


常见问题
------------
如果运行报错，提示缺少依赖（常见于MacOS），这通常是缺少 ``OpenMP`` 库。尝试安装 ``gcc`` 编译器，其中会自带 ``OpenMP``。




