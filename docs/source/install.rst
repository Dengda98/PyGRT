
安装
=============

:Author: Zhu Dengda
:Email:  zhudengda@mail.iggcas.ac.cn

-----------------------------------------------------------


.. |gr| replace:: `Github Releases <https://github.com/Dengda98/PyGRT/releases>`__

依赖
------------

如果想使用 C 程序 :command:`grt` ，

+ `FFTW <https://www.fftw.org/>`_ ，其静态库已链接到预构建版本
+ `NetCDF <https://www.unidata.ucar.edu/software/netcdf>`_ ，其静态库已链接到预构建版本
+ `Seismic Analysis Code (SAC) <http://www.iris.edu/ds/nodes/dmc/forms/sac/>`_ ，需在对应网址申请下载。
  用于用户进一步处理 SAC 格式的输出波形。

如果你想以Python脚本形式使用，依赖已在 :file:`setup.py` 中写好，直接使用 :command:`pip` 安装即可。


安装预构建版本
------------

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

C程序 :command:`grt` 的运行独立于Python，每个C程序可使用 ``-h`` 查看帮助。


从源码构建安装
---------------------

如果安装好预构建版本后，运行 C 或 Python 提示 “GLIBC” 版本缺失以及其它库的版本问题，可尝试从源码从头构建二进制库文件和可执行文件。

1. 安装程序开发所需的基本工具，如 :command:`gcc` 编译器， :command:`make` 工具等。

2. 安装 `NetCDF <https://www.unidata.ucar.edu/software/netcdf>`_ 。

3. 安装 `FFTW <https://www.fftw.org/>`_ 。

  要求安装双精度和单精度两个版本，且要求编译出静态库，并要求在环境变量 :envvar:`LIBRARY_PATH` 中配置静态库路径。如果从源码编译安装 FFTW ，通常在其目录下运行以下命令即可安装+配置路径成功（以 Ubuntu 系统为例）：

  .. code-block:: bash
      
    # 编译双精度版本
    ./configure CFLAGS="-fPIC" 
    make
    make install

    # 编译单精度版本
    ./configure CFLAGS="-fPIC" --enable-float
    make
    make install

4. 构建 **PyGRT** 。

  切换到 **PyGRT** 程序目录（不论是使用 :command:`pip` 安装的还是从 |gr| 下载的），切换到 :rst:dir:`pygrt/C_extension` ，运行 
   
  .. code-block:: bash

      make

  进行构建。成功后会在 :rst:dir:`bin/` 和 :rst:dir:`lib/` 路径下看到新构建出来的可执行文件和库文件。如果正确配置了 :envvar:`PATH` 可尝试运行 :command:`grt -h` 看能否正常打印帮助文档。再运行

  .. code-block:: bash

      make cleanbuild

  可清理构建过程产生的中间文件夹 :rst:dir:`build/` 。

常见问题
------------
+ 如果运行报错，提示缺少依赖（常见于MacOS），这通常是缺少 ``OpenMP`` 库。尝试安装 :command:`gcc` 编译器，其中会自带 ``OpenMP``。
+ “GLIBC” 版本缺失：请尝试从源码构建。




