/**
 * @file   logo.h
 * @author Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
 * @date   2024-12
 * 
 *    logo字符串
 * 
 */

#pragma once


#include <stdio.h>

#include "grt/common/const.h"
#include "grt/common/version.h"
#include "grt/common/colorstr.h"

#if _TEST_WHETHER_WIN32_
#include <windows.h>
#endif


inline GCC_ALWAYS_INLINE void grt_print_logo(){

#if _TEST_WHETHER_WIN32_
    // 在Windows上设置控制台代码页为UTF-8，以显示下方由特殊字符组成的logo
    SetConsoleOutputCP(65001);
#endif

printf(BOLD_GREEN "\n"
"╔═══════════════════════════════════════════════════════════════╗\n"
"║                                                               ║\n"
"║              ██████╗     ██████╗     ████████╗                ║\n"
"║             ██╔════╝     ██╔══██╗    ╚══██╔══╝                ║\n"
"║             ██║  ███╗    ██████╔╝       ██║                   ║\n"
"║             ██║   ██║    ██╔══██╗       ██║                   ║\n"
"║             ╚██████╔╝    ██║  ██║       ██║                   ║\n"
"║              ╚═════╝     ╚═╝  ╚═╝       ╚═╝                   ║\n"
"║                                                               ║\n"
"║                                                               ║\n"
"║               Author: Zhu Dengda                              ║\n"
"║                Email: zhudengda@mail.iggcas.ac.cn             ║\n"
"║        Code Homepage: https://github.com/Dengda98/PyGRT       ║\n"
"║              License: GPL-3.0 license                         ║\n"
"║              Version: %-20s                    ║\n"
"║                                                               ║\n"
"║                                                               ║\n"
"║    A Command-Line Tool for Computing Synthetic Seismograms    ║\n"
"║            in Horizontally Layered Halfspace Model,           ║\n"
"║     using Generalized Reflection-Transmission Method(GRTM)    ║\n"
"║              and Discrete Wavenumber Method(DWM).             ║\n"
"║                                                               ║\n"
"╚═══════════════════════════════════════════════════════════════╝\n" 
DEFAULT_RESTORE, GRT_VERSION);

}