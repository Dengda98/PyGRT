:author: 朱邓达
:date: 2025-04-17

安装
=============

.. |gr| replace:: `Github Releases <https://github.com/Dengda98/PyGRT/releases>`__

依赖
------------

如果想使用 C 程序 :command:`grt` ，

+ |FFTW| *（其静态库已链接到预构建版本）*
+ |NetCDF|  *（其静态库已链接到预构建版本）*
+ `Seismic Analysis Code (SAC) <http://www.iris.edu/ds/nodes/dmc/forms/sac/>`_ ，需在对应网址申请下载。
  用于用户进一步处理 SAC 格式的输出波形。

如果你想以Python脚本形式使用，依赖已在 :file:`setup.py` 中写好，直接使用 :command:`pip` 安装即可。


安装预构建版本
--------------------

+ 目前 **PyGRT** 已分发内置预构建二进制文件的安装包，用户仅需运行以下命令即可（建议使用 `Anaconda <https://anaconda.org>`_ 虚拟环境）

  .. code-block:: bash

      pip install pygrt-kit

  或者从 |gr| 中下载最新版本的程序压缩包（符合自己的操作系统），解压，在根目录下运行以下命令即可

  .. code-block:: bash

      pip install .  

+ 如果你不想使用Python，只想使用传统的命令行形式运行C程序，也可从 |gr| 中下载最新版本的程序压缩包（符合自己的操作系统），
  解压，其中 :rst:dir:`pygrt/C_extension/bin` 和 :rst:dir:`pygrt/C_extension/lib` 为预构建好的可执行文件目录和动态/静态库目录，
  按自己习惯配置环境变量 :envvar:`PATH` 即可（详见下方）。


环境变量配置
-------------
+ 如果你使用 :command:`pip` 安装后，想使用构建好的C程序 :command:`grt` ，需配置环境变量 :envvar:`PATH` 。运行以下命令

  .. code-block:: bash

      python -m pygrt.print

  输出

  .. code-block:: text
    
      PyGRT installation directory: </path/to/installation>
      PyGRT executable file directory: </path/to/installation/bin>
      PyGRT library directory: </path/to/installation/lib>

  将其中的 “PyGRT executable file directory” 路径添加到环境变量 :envvar:`PATH` 中即可。

+ 如果是从 |gr| 上直接下载的压缩包，则只需将解压后的 :rst:dir:`bin/` 路径添加到环境变量 :envvar:`PATH` 中即可。

C程序 :command:`grt` 的运行独立于Python， :command:`grt` 的每个模块可使用 ``-h`` 查看帮助， 例如 :command:`grt greenfn -h` 。


从源码构建安装
---------------------

如有需要，可尝试从源码从头构建二进制库文件和可执行文件。

1. 安装程序开发所需的基本工具，如 :command:`gcc` 编译器， :command:`make` 工具等。

2. 安装 |NetCDF| 。

3. 安装 |FFTW| 。

  要求安装双精度和单精度两个版本的 FFTW，且要求编译出静态库。通常在 FFTW 目录下运行以下命令即可：

  .. code-block:: bash
      
    # 编译双精度版本
    ./configure CFLAGS="-fPIC" 
    make
    sudo make install

    # 编译单精度版本
    ./configure CFLAGS="-fPIC" --enable-float
    make
    sudo make install

  安装好后需确保安装的路径在 :envvar:`LIBRARY_PATH` 中能找到。由于以上运行 configure 时未指定 *--prefix* ，
  默认安装路径一般是 :rst:dir:`/usr/local/` ，因此头文件会在 :rst:dir:`/usr/local/include` 路径下，
  库文件会在 :rst:dir:`/usr/local/lib` 路径下。

4. 构建 **PyGRT** 。

  切换到 **PyGRT** 程序目录（不论是使用 :command:`pip` 安装的还是从 |gr| 下载的），切换到 :rst:dir:`pygrt/C_extension` ，运行 
   
  .. code-block:: bash

      make

  进行构建。对于 Mac 用户，由于 :command:`gcc` 命令对应的不是 GNU 的编译器，因此需要安装 :command:`gcc` 编译器，
  并在以上编译时显式地指定编译器（安装 |FFTW| 的过程也要改），例如

  .. code-block:: bash

      make CC=gcc-14
  
  成功后会在 :rst:dir:`bin/` 和 :rst:dir:`lib/` 路径下看到新构建出来的可执行文件和库文件。如果正确配置了 :envvar:`PATH` 可尝试运行 :command:`grt -h` 看能否正常打印帮助文档。再运行

  .. code-block:: bash

      make cleanbuild

  可清理构建过程产生的中间文件夹 :rst:dir:`build/` 。

常见问题
------------
+ 如果运行报错，提示缺少依赖（常见于MacOS），这通常是缺少 ``OpenMP`` 库。尝试安装 :command:`gcc` 编译器，其中会自带 ``OpenMP``。
+ “GLIBC” 版本缺失：请尝试从源码构建。
+ 从源码构建时提示未找到头文件 ``fftw.h`` : 在环境变量 :envvar:`C_INCLUDE_PATH` 中添加 FFTW 的头文件路径，
  或者在运行 :command:`make` 命令时通过 ``CFLAGS2`` 临时增加 :command:`gcc` 的头文件搜索路径，例如::

    make CFLAGS="-I/usr/local/include -I<其它路径> -I<其它路径>"




