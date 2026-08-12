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
  用于用户进一步处理 SAC 格式的输出波形（可选）。

以 Python 脚本使用时，其余依赖已在 :file:`setup.py` 中写好，直接 :command:`pip` 安装即可。
Python 接口通过调用包内的 :command:`grt` 完成主计算；预构建安装包已按平台内置该可执行文件，
安装后即可使用，**无需再配置** :envvar:`PATH` 。


安装预构建版本
--------------------

目前 **PyGRT** 已在 |gr| 中分发不同平台的预构建安装包（内置对应平台的 :command:`grt` 与库文件）。
用户仅需运行以下命令即可（建议使用 `conda <https://anaconda.org>`_ 虚拟环境）

.. code-block:: bash

    pip install pygrt-kit

或者从 |gr| 中下载符合自己操作系统的程序压缩包，解压后在根目录运行

.. code-block:: bash

    pip install .  

安装完成后即可在 Python 中 ``import pygrt`` 使用。程序会自动定位安装目录内的
:rst:dir:`pygrt/C_extension/bin/grt` ，不必额外配置环境变量。

仅使用命令行 :command:`grt`
--------------------------------
如果你不想使用 Python，只想在终端以命令行形式运行 C 程序，也可从 |gr| 下载对应平台的
``*.tar.gz`` 压缩包（Mac 用户：Apple 芯片选 ``macosx_11_0_arm64`` ，Intel 芯片选 ``macosx_10_9_x86_64`` ）。
解压后，:rst:dir:`pygrt/C_extension/bin` 与 :rst:dir:`pygrt/C_extension/lib` 分别为预构建的可执行文件目录和库目录。
此时需将 :rst:dir:`bin/` 加入环境变量 :envvar:`PATH` ，以便在终端直接调用 :command:`grt` 。

使用 :command:`pip` 安装后若也希望在终端直接运行 :command:`grt` ，可先查看可执行文件路径：

.. code-block:: bash

    python -m pygrt.print

输出形如

.. code-block:: text
  
    PyGRT installation directory: </path/to/installation>
    PyGRT executable file directory: </path/to/installation/bin>
    PyGRT library directory: </path/to/installation/lib>

将其中的 “PyGRT executable file directory” 加入 :envvar:`PATH` 即可。
各模块可用 ``-h`` 查看帮助，例如 :command:`grt greenfn -h` 。


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
  
  成功后会在 :rst:dir:`bin/` 和 :rst:dir:`lib/` 路径下看到新构建出来的可执行文件和库文件。
  Python 侧会自动使用包内刚构建的 :command:`grt` ；若要在终端直接调用，将 :rst:dir:`bin/` 加入 :envvar:`PATH` ，
  并运行 :command:`grt -h` 检查。然后可执行

  .. code-block:: bash

      make cleanbuild

  可清理构建过程产生的中间文件夹 :rst:dir:`build/` 。

常见问题
------------
+ 运行报错，提示缺少依赖（常见于MacOS）
  
  通常是缺少 ``OpenMP`` 库。尝试安装 :command:`gcc` 编译器，其中会自带 ``OpenMP``。

+ “GLIBC” 版本缺失
  
  请尝试从源码构建。

+ 从源码构建时提示未找到头文件 ``fftw.h`` 
  
  在环境变量 :envvar:`C_INCLUDE_PATH` 中添加 FFTW 的头文件路径，
  或者在运行 :command:`make` 命令时通过 ``CFLAGS2`` 临时增加 :command:`gcc` 的头文件搜索路径，例如::

    make CFLAGS="-I/usr/local/include -I<其它路径> -I<其它路径>"
