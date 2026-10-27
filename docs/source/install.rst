:author: 朱邓达
:date: 2025-04-17

安装
=============

.. |gr| replace:: `GitHub Releases <https://github.com/Dengda98/PyGRT/releases>`__

1. 选择安装方式
----------------------------

**PyGRT** 支持在 Python 脚本中调用，也可在终端通过 :command:`grt` 命令运行。
两种方式使用相同的底层计算程序，**命令行方式不需要 Python 环境**。

+ 使用 Python，或同时使用 Python 和命令行：参见 :ref:`install-python`
+ 仅使用命令行：参见 :ref:`install-cli`
+ 修改 C 源码，或预构建版本不适用于当前系统：参见 :ref:`install-source`

预构建版本支持 **Linux（x86_64）、macOS（Intel 和 Apple Silicon）及 Windows（x86_64）**。

.. _install-python:

2. 安装 Python 包
----------------------------

1. 建议在独立环境（如 `conda <https://anaconda.org>`_ 虚拟环境）中安装：

   .. code-block:: bash

       pip install pygrt-kit

   依赖会自动安装。安装时的包名为 **pygrt-kit**，在 Python 中使用时为 **pygrt**。

2. 运行以下命令检查能否正常导入：

   .. code-block:: bash

       python -c "import pygrt; print(pygrt.__version__)"

Python 接口无需配置 :envvar:`PATH`。

2.1 同时使用终端命令（可选）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

无需重复下载，先查询已安装的程序目录：

.. code-block:: bash

    python -m pygrt.print

将输出中 **PyGRT executable file directory** 对应的目录按 :ref:`配置 PATH <install-path>` 加入
:envvar:`PATH`，再运行 :command:`grt -h` 验证。

2.2 从本地安装包安装（可选）
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

也可从 |gr| 下载对应系统和架构的预构建 **.tar.gz** 压缩包，解压后在其根目录运行：

.. code-block:: bash

    pip install .


.. _install-cli:

3. 安装独立命令行程序
----------------------------------------

1. **下载并解压**

   从 |gr| 下载对应系统和架构的预构建 **.tar.gz** 压缩包并解压，保留完整目录。
   可执行文件位于其中的 :file:`pygrt/C_extension/bin/`。

2. **配置 PATH**

   .. _install-path:

   **Linux / macOS**：将下面的路径替换为可执行文件所在目录的绝对路径，在终端运行：

   .. code-block:: bash

       export PATH="/path/to/pygrt/C_extension/bin:$PATH"

   若需长期生效，将这一行加入所用 Shell 的配置文件（例如 Bash 为 :file:`~/.bashrc`，
   Zsh 为 :file:`~/.zshrc`），并重新打开终端。

   **Windows**：在“环境变量”设置中，将可执行文件所在目录的绝对路径添加到用户变量
   **Path**，然后重新打开终端。

3. **验证安装**

   .. code-block:: bash

       grt -h
       grt greenfn -h

   前者列出可用模块，后者查看 :doc:`Module/greenfn` 模块的帮助。
   各模块均按以下格式运行，具体选项见 :doc:`Module/index`：

   .. code-block:: text

       grt <module-name> [<module-options>]


.. _install-source:

4. 从源码构建
----------------------------

仅在需要修改 C 源码或无法使用预构建版本时，才需要自行编译。

1. **获取源码**

   从项目仓库下载源码：

   .. code-block:: bash

       git clone https://github.com/Dengda98/PyGRT.git
       cd PyGRT

2. **安装编译依赖**

   安装支持 OpenMP 的 **GNU GCC**、**make**，以及 |NetCDF| 和 |FFTW| 的开发库
   （包含头文件）。FFTW 需要同时提供双精度和单精度版本，即 **fftw3** 和 **fftw3f**。
   若自行编译这些依赖，需启用 **-fPIC**，以便用于构建共享库。

3. **编译程序**

   在源码根目录下运行：

   .. code-block:: bash

       cd pygrt/C_extension
       make CC=gcc -j

   macOS 自带的 :command:`gcc` 不是 GNU GCC，需另行安装并显式指定名称，
   例如 :command:`make CC=gcc-16 -j`；版本号应与实际安装的编译器一致。

4. **安装与验证**

   根据使用方式完成安装：

   * **命令行**：将生成的 :file:`bin/` 目录按 :ref:`配置 PATH <install-path>` 加入 :envvar:`PATH`，
     运行 :command:`grt -h` 验证
   * **Python**：返回源码根目录，运行 :command:`pip install .`，
     再按 :ref:`install-python` 验证


5. 可选工具
----------------------------------------

部分文档示例使用 `ObsPy <https://docs.obspy.org/>`_ 读取和绘制 SAC 文件的波形，
复现这些示例前请按照对应说明安装 ObsPy。
若习惯使用 `Seismic Analysis Code (SAC) <http://www.iris.edu/ds/nodes/dmc/forms/sac/>`_
处理输出波形，可自行申请下载。运行 **PyGRT** 本身不需要安装 SAC。


6. 常见问题
----------------------------

**终端找不到 grt**
    检查可执行文件所在目录是否已加入 :envvar:`PATH`，配置后重新打开终端。

**Python 提示找不到 grt**
    确认安装的是对应平台的预构建包；若使用仓库源码，需先完成 :ref:`install-source`。

**缺少 OpenMP 运行库（常见于 macOS）**
    安装 GNU GCC 及其 OpenMP 运行库；若仍提示缺少特定版本的库，可在本机从源码构建。

**提示 GLIBC 版本不兼容**
    在当前系统上 :ref:`从源码构建 <install-source>`。

**编译时找不到头文件或库**
    确认已安装 FFTW、NetCDF 的开发库，并将头文件目录加入 :envvar:`C_INCLUDE_PATH`，
    库目录加入 :envvar:`LIBRARY_PATH`。例如，依赖安装在 :file:`/usr/local/` 时：

    .. code-block:: bash

        export C_INCLUDE_PATH="/usr/local/include${C_INCLUDE_PATH:+:$C_INCLUDE_PATH}"
        export LIBRARY_PATH="/usr/local/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"

    然后重新运行编译命令。
