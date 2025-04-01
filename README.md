<p align="center">
  <img src="./figs/logo.png" alt="Image 2" width="45%" />
</p>

<p align="center">
  <img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/Dengda98/PyGRT">
  <img alt="GitHub Actions Workflow Status" src="https://img.shields.io/github/actions/workflow/status/Dengda98/PyGRT/build.yml">
  <img alt="Github Tag" src="https://img.shields.io/github/v/tag/Dengda98/PyGRT">
  <img alt="GitHub License" src="https://img.shields.io/github/license/Dengda98/PyGRT">
</p>

(Detailed documentation is coming soon...)

# Overview
[**PyGRT**](https://github.com/Dengda98/PyGRT): An Efficient and Integrated Python Package for Computing Synthetic Seismograms in a Layered Half-Space Model. 

At present, **PyGRT** can run on
  - [x] Linux
  - [x] macOS
  - [x] Windows

PyGRT is still evolving, and more features will be released in the future.

# Features

- **Dual-Language**:  
  To optimize performance, **PyGRT** uses **C** for its core computational tasks, while **Python** provides a user-friendly interface. Support **script style** and **command line style** to run the program.

- **Parallelization**:  
  Accelerated with **OpenMP** for parallel processing.

- **Integration**:  
  Built on the **Generalized Reflection-Transmission matrix Method (GRTM)** and the **Discrete Wavenumber Method (DWM)**, **PyGRT** integrates the **Peak-Trough Averaging Method (PTAM)** and **Filon’s Integration Method (FIM)** to handle diverse source-receiver distributions. 

- **Modular Design**:   
  Clean and organized code structure, making it easy to extend and maintain.

- **Compatibility**:  
  **PyGRT provides pre-compiled static files**, ensuring ease of installation, usage, and portability across different systems.

<p align="center">
  <img src="./figs/diagram_cut.png" alt="Image 2" width="100%" />
</p>


# Pre-Requisite

+ For Thread-Level Parallel Computing
  - [**OpenMP**](https://www.openmp.org/)
    - For Linux and macOS users: If the GNU compiler is installed on your system, the `OpenMP` library is usually included.  
    - For Windows users: `OpenMP` has been statically linked.  

    **In general, you don't have to worry about it**. However, if the program complains that "`libgomp.so` not found" or "needs more dependencies", you should install `OpenMP`.
    
<br>

+ For Python Script Style  
  - [**Anaconda**](https://anaconda.org) (*recommend*), to build your virtual environment.
  - Other dependencies are declared in `setup.py`, automatically handled by `pip install`.
    
<br>

+ For Command Line Style  
  the output waveforms are binary files in SAC format, you need [**Seismic Analysis Code (SAC)**](http://www.iris.edu/ds/nodes/dmc/forms/sac/) to view and process.
    



# Installation
In **PyGRT**, the C programs and libraries operate independently of Python (not CPython or Cython). If you are not familiar with Python and pip, and prefer the Command Line Style, you can quickly run the program by downloading the latest [**GitHub release**](https://github.com/Dengda98/PyGRT/releases) for your machine. The necessary files are located in the `pygrt/C_extension/bin` and `pygrt/C_extension/lib` folders.

------

Two ways, choose one:
 
1. **PYPI** (recommend)  
  Run the following command in your virtual environment:
    ```Bash
    pip install -v pygrt-kit
    ```


2. [**Github release**](https://github.com/Dengda98/PyGRT/releases)

   - Download the latest [release](https://github.com/Dengda98/PyGRT/releases) for your machine, uncompress, and change the directory.

   - Run the following command in your virtual environment:

      ```bash
      pip install -v .
      ```

3. Build from Source Code.   
   *Not recommend.*

# Setting

For Command Line Style, run
```bash
python -m pygrt.print
```
the outputs are
```
PyGRT installation directory: </path/to/installation>
PyGRT executable file directory: </path/to/installation/bin>
PyGRT library directory: </path/to/installation/lib>
```
and you can 
+ add "executable file directory" to `PATH` environment variable.

Then you can run the command like `grt` in terminal. For each command, use `-h` to see the help message.


# Usage Example

`example/` folder shows some examples in paper. **More examples are coming soon.**

<p align='center'>
  <img alt="multi traces" src="example/multi_traces/multi_traces.png" width="300">
  <img alt="lamb problem" src="example/lamb_problem/stream.png" width="300">
  <img alt="far-field record" src="example/far_field/test.png" width="600">
</p>


# Contact
If you have any questions or suggestions, feel free to reach out:
- **Email**: zhudengda@mail.iggcas.ac.cn
- **GitHub Issues**: You can also raise an issue directly on GitHub.

# Citation

> Zhu D., J. Wang*, J. Hao, S. Yao, Y. Xu, T. Xu and Z. Yao (2025). PyGRT: An Efficient and Integrated Python Package for Computing Synthetic Seismograms in a Layered Half-Space Model. Seismological Research Letters. (under review)
