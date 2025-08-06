Changelog
====================

## PyGRT v0.10.0

This is a **major update** before releasing a _stable version_, primarily concerning the execution style of the C program (the Python interface remains unaffected). Recognizing that compiling multiple executables would **significantly hinder future maintenance and feature expansion**, and considering the current early-stage program isn't too complex yet, I decided to implement this change promptly. 

Inspired by  [`Git`](https://git-scm.com/) and [`GMT`](https://www.generic-mapping-tools.org/) programs, **PyGRT** will now compile into just one executable — **`grt`** — starting with this version. Similar to GMT, it will invoke different modules to perform various computational functions, for example:

``` bash
## Calculate Green's functions (equivalent to the old 'grt' command)
grt greenfn [arg1] [arg2] ...
## Synthesize theoretical seismograms (equivalent to the old 'grt.syn' command)
grt syn [arg1] [arg2] ...
## Calculate static Green's functions (equivalent to the old 'stgrt' command)
grt static greenfn [arg1] [arg2] ...
...
```

**`grt`** and each module's help documentation can be printed using `-h`. The [online documentation](https://pygrt.readthedocs.io/) has been updated accordingly, and all example scripts have been modified to reflect these changes.


* DOC: update installation and remove some unused `*.po` files by [@Dengda98](https://github.com/Dengda98) in [#52](https://github.com/Dengda98/PyGRT/pull/52)
* FIX: allow zero Q^-1, and relevant attenuation scale equals 1.0 by [@Dengda98](https://github.com/Dengda98) in [#53](https://github.com/Dengda98/PyGRT/pull/53)
* Run C-command `grt` like GMT-style by [@Dengda98](https://github.com/Dengda98) in [#55](https://github.com/Dengda98/PyGRT/pull/55)
* FEAT: add grt -v to show the version by [@Dengda98](https://github.com/Dengda98) in [#56](https://github.com/Dengda98/PyGRT/pull/56)
* FIX: update new command in `example/` scripts by [@Dengda98](https://github.com/Dengda98) in [#57](https://github.com/Dengda98/PyGRT/pull/57)
* CI: add CFLAGS="-fPIC" in fftw builds in centos 7 docker by [@Dengda98](https://github.com/Dengda98) in [#58](https://github.com/Dengda98/PyGRT/pull/58)
* STYLE: alyways display `static_<name>` by [@Dengda98](https://github.com/Dengda98) in [#60](https://github.com/Dengda98/PyGRT/pull/60)
* DOC: add changelog.md by [@Dengda98](https://github.com/Dengda98) in [#61](https://github.com/Dengda98/PyGRT/pull/61)
* DOC: change switch color by [@Dengda98](https://github.com/Dengda98) in [#62](https://github.com/Dengda98/PyGRT/pull/62)
* FIX: update new command in `docs/` scripts, and update their description in doc by [@Dengda98](https://github.com/Dengda98) in [#59](https://github.com/Dengda98/PyGRT/pull/59)


**Full Changelog**: [v0.9.2...v0.10.0](https://github.com/Dengda98/PyGRT/compare/v0.9.2...v0.10.0)


##   PyGRT v0.9.2

In Github Actions, use CentOS 7 docker to compile PyGRT, to support `glibc` >= 2.17 (ensure compatibility as much as possible).

* REFAC: some tiny change for downward compatibility by [@Dengda98](https://github.com/Dengda98) in [#50]([#50](https://github.com/Dengda98/PyGRT/pull/50))
* Bump to PyGRT v0.9.2 by [@Dengda98](https://github.com/Dengda98) in [#51](https://github.com/Dengda98/PyGRT/pull/51)


**Full Changelog**: [v0.9.1...v0.9.2](https://github.com/Dengda98/PyGRT/compare/v0.9.1...v0.9.2)


##   PyGRT v0.9.1

Bugfix release addressing issues in **v0.9.0**'s source located in liquid. Update recommended: `pip install --upgrade pygrt-kit`

* DOC: add note about near-field in FIM and SAFIM by [@Dengda98](https://github.com/Dengda98) in [#46](https://github.com/Dengda98/PyGRT/pull/46)
* FEAT: update colorstr by [@Dengda98](https://github.com/Dengda98) in [#47](https://github.com/Dengda98/PyGRT/pull/47)
* FIX: add some warning about source in liquid by [@Dengda98](https://github.com/Dengda98) in [#48](https://github.com/Dengda98/PyGRT/pull/48)
* Bump to PyGRT v0.9.1 by [@Dengda98](https://github.com/Dengda98) in [#49](https://github.com/Dengda98/PyGRT/pull/49)


**Full Changelog**: [v0.9.0...v0.9.1](https://github.com/Dengda98/PyGRT/compare/v0.9.0...v0.9.1)


##  PyGRT v0.9.0

I'm very excited to announce the new release of **PyGRT v0.9.0** has supported for computing synthetic seismograms with **liquid layers**, with results validated against CPS330's `rspec96` module for accuracy. The [documentation](https://pygrt.readthedocs.io/zh-cn/dev/Formula/index.html) has been updated to include the corresponding formulas. This release is particularly useful for studies involving ocean. 

Upgrade via `pip install --upgrade pygrt-kit`.



*  by [@Dengda98](https://github.com/Dengda98) in [#29](https://github.com/Dengda98/PyGRT/pull/29)
* BUILD: fix makefile by [@Dengda98](https://github.com/Dengda98) in [#30](https://github.com/Dengda98/PyGRT/pull/30)
* REFAC: update temp var for `R_EVL` to improve code readability by [@Dengda98](https://github.com/Dengda98) in [#31](https://github.com/Dengda98/PyGRT/pull/31)
* DOC: add RT.rst and RT_formula.ipynb by [@Dengda98](https://github.com/Dengda98) in [#33](https://github.com/Dengda98/PyGRT/pull/33)
* Support liquid layer by [@Dengda98](https://github.com/Dengda98) in [#34](https://github.com/Dengda98/PyGRT/pull/34)
* add document about R/T matrix formula by [@Dengda98](https://github.com/Dengda98) in [#35](https://github.com/Dengda98/PyGRT/pull/35)
* TEST: add compare_results_liquid example by [@Dengda98](https://github.com/Dengda98) in [#36](https://github.com/Dengda98/PyGRT/pull/36)
* FIX: use nonzero minimum Vb in python by [@Dengda98](https://github.com/Dengda98) in [#37](https://github.com/Dengda98/PyGRT/pull/37)
* allow zero vs for computing sythetic seismogram and stress  by [@Dengda98](https://github.com/Dengda98) in [#38](https://github.com/Dengda98/PyGRT/pull/38)
* FIX: raise error if receiver located on the interface where liquid exist by [@Dengda98](https://github.com/Dengda98) in [#39](https://github.com/Dengda98/PyGRT/pull/39)
* TEST: add seafloor example by [@Dengda98](https://github.com/Dengda98) in [#40](https://github.com/Dengda98/PyGRT/pull/40)
* DOC: add a note to support liquid layer by [@Dengda98](https://github.com/Dengda98) in [#41](https://github.com/Dengda98/PyGRT/pull/41)
* REFAC: remove unused code by [@Dengda98](https://github.com/Dengda98) in [#42](https://github.com/Dengda98/PyGRT/pull/42)
* DOC: add description about liquid layer supporting by [@Dengda98](https://github.com/Dengda98) in [#43](https://github.com/Dengda98/PyGRT/pull/43)
* Bump to PyGRT v0.9.0 by [@Dengda98](https://github.com/Dengda98) in [#44](https://github.com/Dengda98/PyGRT/pull/44)
* CI: add some tests by [@Dengda98](https://github.com/Dengda98) in [#45](https://github.com/Dengda98/PyGRT/pull/45)


**Full Changelog**: [v0.8.0...v0.9.0](https://github.com/Dengda98/PyGRT/compare/v0.8.0...v0.9.0)


##  PyGRT v0.8.0

I'm excited to announce that **now PyGRT supports Self-Adaptive Filon's Integration Method(SAFIM)**, a powerful tool to efficiently compute Green's Functions at large epicentral distances.

![image](https://pygrt.readthedocs.io/zh-cn/v0.8.0/_images/safim.png)

* Refactor codes, and support Self-Adaptive Filon's Integration Method (SAFIM) by [@Dengda98](https://github.com/Dengda98) in [#27](https://github.com/Dengda98/PyGRT/pull/27)
* DOC: add SAFIM by [@Dengda98](https://github.com/Dengda98) in [#28](https://github.com/Dengda98/PyGRT/pull/28)


**Full Changelog**: [v0.7.0...v0.8.0](https://github.com/Dengda98/PyGRT/compare/v0.7.0...v0.8.0)


##  PyGRT v0.7.0


* add English docs by [@Dengda98](https://github.com/Dengda98) in [#24](https://github.com/Dengda98/PyGRT/pull/24)
* REMOVE: remove `iwk0` parameter by [@Dengda98](https://github.com/Dengda98) in [#25](https://github.com/Dengda98/PyGRT/pull/25)
* update FIM (code and doc), add warning about `dk` in FIM by [@Dengda98](https://github.com/Dengda98) in [#26](https://github.com/Dengda98/PyGRT/pull/26)


**Full Changelog**: [v0.6.0...v0.7.0](https://github.com/Dengda98/PyGRT/compare/v0.6.0...v0.7.0)


##  PyGRT v0.6.0

* DOC: add "strain and stress calculation" in Tutorial, and fix some typos by [@Dengda98](https://github.com/Dengda98) in [#16](https://github.com/Dengda98/PyGRT/pull/16)
* Update stats file (kernel functions value in integration ) format, python reading functions, add `grt.k2a` command by [@Dengda98](https://github.com/Dengda98) [#17](https://github.com/Dengda98/PyGRT/pull/17)
* DOC: add `integ_converg` in Toturial by [@Dengda98](https://github.com/Dengda98) [#18](https://github.com/Dengda98/PyGRT/pull/18)
* FEAT: add dist in output of function `pygrt.utils.read_statsfile_ptam` by [@Dengda98](https://github.com/Dengda98) [#19](https://github.com/Dengda98/PyGRT/pull/19)
* DOC: Refactor tutorial scripts and documentation for clarity and consistency by [@Dengda98](https://github.com/Dengda98) [#20](https://github.com/Dengda98/PyGRT/pull/20)
* DOC: add `grt.travt` in Tutorial by [@Dengda98](https://github.com/Dengda98) [#21](https://github.com/Dengda98/PyGRT/pull/21)
* add `read_kernels_freqs()` function, add "kernel frequency response" in Tutorial by [@Dengda98](https://github.com/Dengda98) [#22](https://github.com/Dengda98/PyGRT/pull/22)
* Support rotation tensor calculation in both dynamic and static case by [@Dengda98](https://github.com/Dengda98) [#23](https://github.com/Dengda98/PyGRT/pull/23)

**Full Changelog**: [v0.5.1...v0.6.0](https://github.com/Dengda98/PyGRT/compare/v0.5.1...v0.6.0)

![view_stats](https://github.com/Dengda98/PyGRT/raw/main/example/view_integ_stats/view_stats.png)
![imag_G](https://github.com/Dengda98/PyGRT/raw/main/example/kernel_freq_response/imag_G.png)


##  PyGRT v0.5.1

+ **Now both Python and C can compute displacement, strain and stress in dynamic and static cases.**
+ **Fix some math error and text error**, such as, unit mixing (`km` and `cm`) in spatial derivatives in Python, coordinate of moment tensor is NED (North, East, Downward)  , ...
+ **Chinese document is available**, see [here](https://pygrt.readthedocs.io/zh-cn/). Thanks @xichaoqiang for revision.


* Support static displacements, strain and stress calculations in Python by [@Dengda98](https://github.com/Dengda98) in [#14](https://github.com/Dengda98/PyGRT/pull/14)
* Build online docs, and upload to ReadtheDocs by [@Dengda98](https://github.com/Dengda98) in [#15](https://github.com/Dengda98/PyGRT/pull/15)
> `docs/` folder has been removed from the following `*.tar.gz`

**Full Changelog**: [v0.4.0...v0.5.1](https://github.com/Dengda98/PyGRT/compare/v0.4.0...v0.5.1)



## PyGRT v0.4.0

Now in C-level, **PyGRT** can compute **displacements, strain and stress in both dynamic and static case.** Python-level will follow up in later minor version.

* Rearrange directory of C source files by [@Dengda98](https://github.com/Dengda98) in [#6](https://github.com/Dengda98/PyGRT/pull/6)
* FIX: segfault when call `write_stats` in `ptam.c` by [@Dengda98](https://github.com/Dengda98) in [#7](https://github.com/Dengda98/PyGRT/pull/7)
* FEAT: support spatial derivatives of displacements by [@Dengda98](https://github.com/Dengda98) in [#9](https://github.com/Dengda98/PyGRT/pull/9)
* add simple tool `grt.b2a` to convert waveforms data in SAC file into ASCII file. by [@Dengda98](https://github.com/Dengda98) in [#10](https://github.com/Dengda98/PyGRT/pull/10)
* Add some features and fix some bugs about source signals by [@Dengda98](https://github.com/Dengda98) in [#11](https://github.com/Dengda98/PyGRT/pull/11)
* Support stress and strain calculation by [@Dengda98](https://github.com/Dengda98) in [#12](https://github.com/Dengda98/PyGRT/pull/12)
* Big update: Support static displacements calculation, Support strain and stress calculation for dynamic and static case by [@Dengda98](https://github.com/Dengda98) in [#13](https://github.com/Dengda98/PyGRT/pull/13)


**Full Changelog**: [v0.2.0...v0.4.0](https://github.com/Dengda98/PyGRT/compare/v0.2.0...v0.4.0)


## PyGRT v0.2.0 

Support Linux, MacOS and Windows

- [x] Linux
- [x] MacOS 
- [x] **Windows**

* Pre-compiled binary programs and Libraries on Linux, MacOS and **Windows**. You can simply run `pip install pygrt-kit` to install the new version.
  - Linux and MacOS. 
    Dynamic library of `OpenMP` is needed, while in general it has already been included in `GNU` compiler. So if program complain that "libgomp.so not found" or python said "you may need other dependency", just add one step to install `OpenMP`.
  - **Native Windows!!**  
    Different from Linux and MacOS, `OpenMP` are statically linked. You can just open a `cmd` (not `PowerShell`) or `Anaconda Prompt`, directly run the program.

* Improve compatibility, and fix some bugs by [@Dengda98](https://github.com/Dengda98) in [#4](https://github.com/Dengda98/PyGRT/pull/4)



**Full Changelog**: [v0.2.0](https://github.com/Dengda98/PyGRT/commits/v0.2.0)