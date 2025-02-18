<p align="center">
  <img src="./figs/logo.png" alt="Image 2" width="45%" />
</p>

<p align="center">
  <img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/Dengda98/PyGRT">
  <img alt="GitHub License" src="https://img.shields.io/github/license/Dengda98/PyGRT">
</p>



(Detailed documentation is coming soon...)

# Overview
**PyGRT**: An Efficient and Integrated Python Package for Computing Synthetic Seismograms in a Layered Half-Space Model. 

# Features

- **Dual-Language**:  
  To optimize performance, **PyGRT** uses **C** for its core computational tasks, while **Python** provides a user-friendly interface. Support **script style** and **command line style** to run the program.

- **Parallelization**:  
  Accelerated with **OpenMP** for parallel processing.

- **Integration**:  
  Built on the **Generalized Reflection-Transmission matrix Method (GRTM)** and the **Discrete Wavenumber Method (DWM)**, **PyGRT** integrates the **Peak-Trough Averaging Method (PTAM)** and **Filon’s Integration Method (FIM)** to handle diverse source-receiver distributions. 

- **Modular Design**:   
  Clean and organized code structure, making it easy to extend and maintain.


<p align="center">
  <img src="./figs/diagram_cut.png" alt="Image 2" width="100%" />
</p>


# pre-requisite

1. [Anaconda](https://anaconda.org) for python script style in virtual environment.
2. [FFTW](https://fftw.org/) for command line style.

# Installation

Two ways, choose one:
 
1. **PYPI** (recommend)  
  Run the following command in your virtual environment:
   ```Bash
   pip install -v pygrt-kit
   ```


2. Github 

   - Download the latest [release](https://github.com/Dengda98/PyGRT/releases), uncompress, and change the diretory.

   - Run the following command in your virtual environment:

      ```bash
      pip install -v .
      ```

If you're more comfortable with the command line style, run
```bash
python -m pygrt.print
```
the outputs are
```
PyGRT installation directory: ...
PyGRT executable file directory: ...
PyGRT dynamic library directory: ...
```
and you can 
+ add "executable file directory" it to `PATH` variable.
+ add "dynamic library directory" it to `LD_LIBRARY_PATH` variable.

Then you can run the command like `grt` in terminal. For each command, use `-h` to see the help message.


# Usage Example

`example/` folder shows some examples in paper. **More examples are coming soon.**


# Contact
If you have any questions or suggestions, feel free to reach out:
- **Email**: zhudengda@mail.iggcas.ac.cn
- **GitHub Issues**: You can also raise an issue directly on GitHub.

# Citation

> Zhu D., J. Wang*, J. Hao, S. Yao, Y. Xu, T. Xu and Z. Yao (2025). PyGRT: An Efficient and Integrated Python Package for Computing Synthetic Seismograms in a Layered Half-Space Model. Seismological Research Letters. (submitted)